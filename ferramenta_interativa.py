"""
FERRAMENTA INTERATIVA DE ANÁLISE SCG
=====================================
Usa aceleração média (AccMag) + limpeza de ruído.
Interface interativa para:
  JANELA 1: Marcar pontos fiduciais na onda média SCG
  JANELA 2: Navegar e validar ciclos individuais
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Backend interativo (Tkinter)
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons
from scipy import signal
import xlrd
from scipy.interpolate import interp1d
import xlrd
from pathlib import Path
import pickle
import warnings
warnings.filterwarnings('ignore')

plt.rcParams['figure.dpi'] = 100
plt.rcParams['toolbar'] = 'none'

DATA_DIR = Path(r"C:\Users\kelse\OneDrive\Área de Trabalho\DADOS_PIBIC")
SCG_DIR = DATA_DIR / "PESQUISA_PIBIC_SISMO"
FS = 100.0

# ============================================================
# PONTOS FIDUCIAIS (Sørensen et al., 2018, Sci. Rep. 8:15455)
# 8 pontos com correlação fisiológica validada por ultrassom.
# ============================================================
# nome -> (evento fisiológico, correlação r, cor, complexo, alias válvula)
FIDUCIAL_DEFS = {
    'Bs': ('Sístole atrial (início onda A)', 0.75, '#e6194B', 's', None),
    'Cs': ('Pico de influxo atrial',          0.63, '#f58231', 's', None),
    'Es': ('Fechamento mitral',               0.71, '#ffbe0b', 's', 'MC'),
    'Gs': ('Abertura aórtica',                0.60, '#3cb44b', 's', 'AO'),
    'Ks': ('Pico de ejeção sistólica',        0.42, '#42d4f4', 's', None),
    'Bd': ('Fechamento aórtico',              0.94, '#4363d8', 'd', 'AC'),
    'Fd': ('Abertura mitral',                 0.87, '#911eb4', 'd', 'MO'),
    'Gd': ('Pico onda E (enchimento vent.)',  0.79, '#f032e6', 'd', None),
}
# Mapa alias-de-válvula -> ponto fiducial (compatibilidade retroativa)
ALIAS_PARA_FIDUCIAL = {v[4]: k for k, v in FIDUCIAL_DEFS.items() if v[4]}

# ============================================================
# FUNÇÕES DE PROCESSAMENTO
# ============================================================

def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = fs / 2
    b, a = signal.butter(order, [lowcut / nyq, highcut / nyq], btype='band')
    return b, a

def apply_filter(data, b, a):
    pad = min(len(data) // 10, 500)
    padded = np.pad(data, pad, mode='reflect')
    filtered = signal.filtfilt(b, a, padded)
    return filtered[pad:-pad]

def remove_noise_sqi(sinal, fs):
    """
    Limpeza de ruído usando filtragem + remoção de artefatos.
    Pipeline:
    1. Filtro passa-banda 2-39 Hz
    2. Remove outliers de amplitude (>4 desvios padrão)
    3. Suavização adaptativa
    """
    b, a = butter_bandpass(2, 39, fs, 4)
    sinal_filt = apply_filter(sinal, b, a)
    
    # Remover regiões com amplitude anômala
    threshold = 4 * np.std(sinal_filt)
    mask_bad = np.abs(sinal_filt) > threshold
    
    # Interpolar pontos ruins (manter como NaN para ignorar)
    sinal_clean = sinal_filt.copy()
    sinal_clean[mask_bad] = np.nan
    
    # Interpolação linear para preencher NaNs
    if np.any(np.isnan(sinal_clean)):
        x = np.arange(len(sinal_clean))
        good = ~np.isnan(sinal_clean)
        if np.sum(good) > 2:
            sinal_clean = np.interp(x, x[good], sinal_clean[good])
        else:
            sinal_clean = np.nan_to_num(sinal_clean)
    
    return sinal_clean

def estimate_fc_autocorr(sinal, fs):
    """Estima FC via autocorrelação do envelope."""
    analytic = signal.hilbert(sinal)
    env = np.abs(analytic)
    nyq = fs / 2
    b, a = signal.butter(4, 5 / nyq, btype='low')
    env_filt = signal.filtfilt(b, a, env)
    
    n = len(env_filt)
    acorr = signal.correlate(env_filt - np.mean(env_filt), env_filt - np.mean(env_filt), mode='full', method='auto')
    acorr = acorr / acorr[n-1]
    half = acorr[n-1:]
    half_lags = np.arange(len(half)) / fs
    
    mask = (half_lags >= 0.35) & (half_lags <= 1.5)
    if not np.any(mask):
        return np.nan
    
    peaks, props = signal.find_peaks(half, height=0.05, distance=int(0.15*fs), prominence=0.03)
    valid = peaks[(half_lags[peaks] >= 0.35) & (half_lags[peaks] <= 1.5)]
    
    if len(valid) == 0:
        return np.nan
    
    best = valid[np.argmax(half[valid])]
    return 60.0 / half_lags[best]

def detect_initial_beats(sinal_clean, fs):
    """Detecção preliminar de batimentos (conservadora)."""
    fc = estimate_fc_autocorr(sinal_clean, fs)
    if np.isnan(fc):
        fc = 75.0
    
    min_dist = max(0.35, 60.0 / fc * 0.7)
    threshold = np.std(sinal_clean) * 0.5
    
    peaks, props = signal.find_peaks(
        sinal_clean,
        height=threshold,
        distance=int(min_dist * fs),
        prominence=threshold * 0.3
    )
    
    return peaks, fc

def extract_cycles(sinal, peaks, fs, before=0.3, after=0.6):
    """Extrai ciclos SCG alinhados pelos picos."""
    n_before = int(before * fs)
    n_after = int(after * fs)
    cycles = []
    valid_peaks = []
    
    for p in peaks:
        if p - n_before >= 0 and p + n_after < len(sinal):
            cycles.append(sinal[p - n_before:p + n_after])
            valid_peaks.append(p)
    
    return cycles, valid_peaks

# ============================================================
# FERRAMENTA INTERATIVA
# ============================================================

class SCGInteractiveTool:
    """Ferramenta interativa para análise SCG."""
    
    def __init__(self, sinal, t, cycles, peaks, patient_num, grupo, idade):
        self.sinal = sinal
        self.t = t
        self.cycles = cycles
        self.peaks = peaks
        self.patient_num = patient_num
        self.grupo = grupo
        self.idade = idade
        
        # Comprimento do ciclo
        if len(cycles) > 0:
            self.cycle_len = len(cycles[0])
            self.time_axis = np.linspace(-0.3, 0.6, self.cycle_len)
        else:
            self.cycle_len = 90
            self.time_axis = np.linspace(-0.3, 0.6, 90)
        
        # Estado
        self.fiducial_points = {}  # {nome: índice no template}
        self.cycle_status = {}  # {idx: 'good'/'bad'/'unmarked'}
        self.current_cycle = 0
        self.corrected_peaks = {}  # {idx: índice no sinal}
        
        # Onda média
        if len(cycles) > 2:
            self.ensemble = np.mean(self.cycles, axis=0)
            self.ensemble_std = np.std(self.cycles, axis=0)
        else:
            self.ensemble = np.zeros(self.cycle_len)
            self.ensemble_std = np.zeros(self.cycle_len)
        
        self.setup_ui()
    
    def setup_ui(self):
        """Cria as duas janelas."""
        # ========== JANELA 1: Onda Média ==========
        self.fig1, (self.ax1, self.ax1_info) = plt.subplots(2, 1, figsize=(14, 10),
                                                              gridspec_kw={'height_ratios': [4, 1]})
        self.fig1.canvas.manager.set_window_title(f'Paciente #{self.patient_num} - Onda Média SCG')
        self.fig1.subplots_adjust(right=0.80)
        
        # Plot onda média
        self.ax1.plot(self.time_axis, self.ensemble, 'b-', linewidth=2.5, label='Onda Média')
        self.ax1.fill_between(self.time_axis, self.ensemble - self.ensemble_std,
                              self.ensemble + self.ensemble_std, color='blue', alpha=0.15)
        self.ax1.axvline(x=0, color='red', linestyle='--', alpha=0.5, linewidth=1.5)
        self.ax1.set_title(f'Paciente #{self.patient_num} | {self.grupo} | {self.idade} anos | n={len(self.cycles)} ciclos',
                           fontsize=14, fontweight='bold')
        self.ax1.set_xlabel('Tempo (s)', fontsize=11)
        self.ax1.set_ylabel('Aceleração (m/s²)', fontsize=11)
        self.ax1.grid(True, alpha=0.3)
        self.ax1.legend(loc='upper right')
        
        # Pontos fiduciais já marcados (8 pontos do artigo)
        self.fiducial_scatter = {}
        self.fiducial_labels = {}
        for nome, (_evento, _r, cor, _cx, alias) in FIDUCIAL_DEFS.items():
            leg = f'{nome} ({alias})' if alias else nome
            self.fiducial_scatter[nome] = self.ax1.scatter([], [], c=cor, s=130,
                                                            marker='o', zorder=5,
                                                            edgecolors='black', linewidths=0.7,
                                                            label=leg)
            self.fiducial_labels[nome] = self.ax1.annotate(
                '', xy=(0, 0), xytext=(0, 10), textcoords='offset points',
                ha='center', fontsize=9, fontweight='bold', color=cor,
                zorder=6, visible=False)
        self.ax1.legend(loc='upper right', fontsize=8, ncol=2)
        
        # Painel de texto
        self.ax1_info.axis('off')
        self.info_text = self.ax1_info.text(0.02, 0.5, self._get_info_text(),
                                            fontsize=10, family='monospace', va='center')
        
        # Seletor de ponto fiducial (RadioButtons) na lateral direita
        self.radio_labels = []
        for nome, (_evento, _r, _cor, _cx, alias) in FIDUCIAL_DEFS.items():
            self.radio_labels.append(f'{nome} · {alias}' if alias else f'{nome}')
        ax_radio = plt.axes([0.82, 0.30, 0.17, 0.55])
        ax_radio.set_title('Ponto a marcar', fontsize=10, fontweight='bold')
        self.radio = RadioButtons(ax_radio, self.radio_labels, active=0)
        for lbl, (nome, defs) in zip(self.radio.labels, FIDUCIAL_DEFS.items()):
            lbl.set_color(defs[2])
            lbl.set_fontsize(9)
        self.radio.on_clicked(self._on_radio_select)
        # Ponto atualmente selecionado
        self.current_fiducial = list(FIDUCIAL_DEFS.keys())[0]

        # Botão SALVAR
        ax_save = plt.axes([0.82, 0.20, 0.17, 0.05])
        self.btn_save = Button(ax_save, 'SALVAR', color='lightgreen')
        self.btn_save.on_clicked(self._save_and_next)

        self.fig1.canvas.mpl_connect('button_press_event', self._on_click_fiducial)
        
        # ========== JANELA 2: Ciclos Individuais ==========
        self.fig2, (self.ax2, self.ax2_ctrl) = plt.subplots(2, 1, figsize=(14, 10),
                                                              gridspec_kw={'height_ratios': [4, 1]})
        self.fig2.canvas.manager.set_window_title(f'Paciente #{self.patient_num} - Ciclos Individuais')
        
        self.ax2.set_title('Ciclo 0', fontsize=13)
        self.ax2.set_xlabel('Tempo (s)')
        self.ax2.set_ylabel('Aceleração (m/s²)')
        self.ax2.grid(True, alpha=0.3)
        
        # Linha do ciclo atual
        self.cycle_line, = self.ax2.plot([], [], 'b-', linewidth=2, label='Ciclo atual')
        self.cycle_ensemble, = self.ax2.plot(self.time_axis, self.ensemble, 'gray', 
                                             linewidth=1, alpha=0.5, label='Onda média')
        self.ax2.axvline(x=0, color='red', linestyle='--', alpha=0.3)
        self.ax2.legend()
        
        # Marcação de pico
        self.peak_marker = None
        self.marking_peak = False
        
        # Status
        self.ax2_ctrl.axis('off')
        self.status_text = self.ax2_ctrl.text(0.02, 0.5, '', fontsize=11, 
                                               family='monospace', va='center')
        
        # Botões Janela 2
        ax_prev = plt.axes([0.15, 0.02, 0.10, 0.04])
        ax_next = plt.axes([0.27, 0.02, 0.10, 0.04])
        ax_good = plt.axes([0.40, 0.02, 0.10, 0.04])
        ax_bad = plt.axes([0.52, 0.02, 0.10, 0.04])
        ax_mark_peak = plt.axes([0.65, 0.02, 0.12, 0.04])
        ax_done = plt.axes([0.80, 0.02, 0.10, 0.04])
        
        self.btn_prev = Button(ax_prev, '<< Anterior')
        self.btn_next = Button(ax_next, 'Próximo >>')
        self.btn_good = Button(ax_good, '✓ Bom', color='lightgreen')
        self.btn_bad = Button(ax_bad, '✗ Ruim', color='lightcoral')
        self.btn_mark_peak = Button(ax_mark_peak, 'Marcar Pico')
        self.btn_done = Button(ax_done, 'CONCLUIR', color='gold')
        
        self.btn_prev.on_clicked(self._prev_cycle)
        self.btn_next.on_clicked(self._next_cycle)
        self.btn_good.on_clicked(lambda e: self._classify_cycle('good'))
        self.btn_bad.on_clicked(lambda e: self._classify_cycle('bad'))
        self.btn_mark_peak.on_clicked(self._toggle_peak_marking)
        self.btn_done.on_clicked(self._finish)
        
        self.fig2.canvas.mpl_connect('button_press_event', self._on_click_cycle)
        
        self._update_cycle_display()
        self._update_status()
    
    def _get_info_text(self):
        text = f"Paciente #{self.patient_num} | Grupo: {self.grupo} | Idade: {self.idade}\n"
        text += f"Ciclos detectados: {len(self.cycles)}\n\n"
        text += "INSTRUÇÕES:\n"
        text += "1. Selecione o ponto na lista à direita e clique no gráfico\n"
        text += "2. 8 pontos fiduciais (Sørensen 2018): Bs Cs Es Gs Ks | Bd Fd Gd\n"
        text += "3. Clique SALVAR quando terminar\n\n"
        text += "PONTOS MARCADOS:\n"
        if self.fiducial_points:
            for nome in FIDUCIAL_DEFS:
                if nome in self.fiducial_points:
                    idx = self.fiducial_points[nome]
                    t_fid = self.time_axis[idx]
                    alias = FIDUCIAL_DEFS[nome][4]
                    tag = f" ({alias})" if alias else ""
                    text += f"  {nome}{tag}: t={t_fid:+.3f}s\n"
        else:
            text += "  (nenhum)\n"
        return text

    def _on_radio_select(self, label):
        """Seleciona o ponto fiducial a partir do rótulo do RadioButton."""
        nome = label.split(' ')[0]
        self.current_fiducial = nome
        evento, r, _cor, cx, alias = FIDUCIAL_DEFS[nome]
        alias_txt = f" (={alias})" if alias else ""
        self.info_text.set_text(
            self._get_info_text()
            + f"\n>>> Clique no gráfico para marcar {nome}{alias_txt}\n"
            + f"    {evento} | complexo {cx} | r={r:.2f} <<<")
        self.fig1.canvas.draw_idle()

    def _start_marking(self, nome):
        self.current_fiducial = nome
        self.info_text.set_text(self._get_info_text() + f"\n>>> Clique no gráfico para marcar {nome} <<<")
        self.fig1.canvas.draw_idle()
    
    def _on_click_fiducial(self, event):
        if self.current_fiducial is None:
            return
        if event.inaxes != self.ax1:
            return
        
        # Encontrar índice mais próximo
        t_click = event.xdata
        idx = np.argmin(np.abs(self.time_axis - t_click))
        
        nome = self.current_fiducial
        self.fiducial_points[nome] = idx
        
        # Atualizar scatter
        x, y = self.time_axis[idx], self.ensemble[idx]
        self.fiducial_scatter[nome].set_offsets([[x, y]])
        # Atualizar rótulo de texto sobre o ponto
        lbl = self.fiducial_labels[nome]
        lbl.set_text(nome)
        lbl.xy = (x, y)
        lbl.set_visible(True)
        
        self.info_text.set_text(self._get_info_text())
        self.fig1.canvas.draw_idle()
    
    def _update_cycle_display(self):
        if len(self.cycles) == 0:
            return
        
        ciclo = self.cycles[self.current_cycle]
        self.cycle_line.set_data(self.time_axis, ciclo)
        
        # Ajustar limites y
        y_min = min(np.min(ciclo), np.min(self.ensemble)) - 0.02
        y_max = max(np.max(ciclo), np.max(self.ensemble)) + 0.02
        self.ax2.set_ylim(y_min, y_max)
        
        status = self.cycle_status.get(self.current_cycle, '---')
        self.ax2.set_title(f'Ciclo {self.current_cycle+1}/{len(self.cycles)}  [{status}]', 
                          fontsize=13, fontweight='bold')
        
        self.fig2.canvas.draw_idle()
    
    def _prev_cycle(self, event):
        if self.current_cycle > 0:
            self.current_cycle -= 1
            self._update_cycle_display()
            self._update_status()
    
    def _next_cycle(self, event):
        if self.current_cycle < len(self.cycles) - 1:
            self.current_cycle += 1
            self._update_cycle_display()
            self._update_status()
    
    def _classify_cycle(self, status):
        self.cycle_status[self.current_cycle] = status
        self._update_cycle_display()
        self._update_status()
        # Auto-avançar
        if self.current_cycle < len(self.cycles) - 1:
            self.current_cycle += 1
            self._update_cycle_display()
            self._update_status()
    
    def _toggle_peak_marking(self, event):
        self.marking_peak = not self.marking_peak
        if self.marking_peak:
            self.btn_mark_peak.color = 'yellow'
            self.btn_mark_peak.label.set_text('CLIQUE no pico...')
            self.status_text.set_text(">>> Clique no gráfico onde está o pico AO real <<<")
        else:
            self.btn_mark_peak.color = 'lightgray'
            self.btn_mark_peak.label.set_text('Marcar Pico')
        self.fig2.canvas.draw_idle()
    
    def _on_click_cycle(self, event):
        if not self.marking_peak:
            return
        if event.inaxes != self.ax2:
            return
        
        t_click = event.xdata
        idx_in_cycle = np.argmin(np.abs(self.time_axis - t_click))
        
        # Converter para índice no sinal original
        peak_idx = self.peaks[self.current_cycle]
        cycle_start = peak_idx - self.cycle_len // 2
        corrected_idx = cycle_start + idx_in_cycle
        
        self.corrected_peaks[self.current_cycle] = corrected_idx
        
        # Mostrar marcador
        if self.peak_marker:
            self.peak_marker.remove()
        self.peak_marker = self.ax2.scatter([t_click], [self.cycles[self.current_cycle][idx_in_cycle]],
                                            c='red', s=150, marker='*', zorder=10)
        
        self.marking_peak = False
        self.btn_mark_peak.color = 'lightgray'
        self.btn_mark_peak.label.set_text('Marcar Pico')
        self._update_status()
        self.fig2.canvas.draw_idle()
    
    def _update_status(self):
        good = sum(1 for v in self.cycle_status.values() if v == 'good')
        bad = sum(1 for v in self.cycle_status.values() if v == 'bad')
        unmarked = len(self.cycles) - good - bad
        corrected = len(self.corrected_peaks)
        
        text = f"Ciclos: {good} bons | {bad} ruins | {unmarked} não avaliados | {corrected} picos corrigidos\n"
        text += f"Atual: Ciclo {self.current_cycle+1}/{len(self.cycles)}\n"
        self.status_text.set_text(text)
        self.fig2.canvas.draw_idle()
    
    def _save_and_next(self, event):
        """Salva anotações e fecha as janelas."""
        self.save_annotations()
        print(f"\n✅ Anotações salvas para paciente #{self.patient_num}!")
        print(f"   Pontos fiduciais: {self.fiducial_points}")
        print(f"   Ciclos bons: {sum(1 for v in self.cycle_status.values() if v=='good')}")
        print(f"   Picos corrigidos: {len(self.corrected_peaks)}")
        plt.close('all')
    
    def _finish(self, event):
        self._save_and_next(event)
    
    def save_annotations(self):
        out_path = DATA_DIR / f"anotacoes_paciente_{self.patient_num}.pkl"
        # Aliases de válvula (AO/AC/MO/MC) para compatibilidade retroativa
        alias_points = {}
        alias_times = {}
        for alias, nome in ALIAS_PARA_FIDUCIAL.items():
            if nome in self.fiducial_points:
                idx = self.fiducial_points[nome]
                alias_points[alias] = idx
                alias_times[alias] = self.time_axis[idx]
        data = {
            'patient': self.patient_num,
            'grupo': self.grupo,
            'idade': self.idade,
            'fiducial_points': self.fiducial_points,
            'fiducial_times': {k: self.time_axis[v] for k, v in self.fiducial_points.items()},
            'fiducial_defs': {k: {'evento': v[0], 'r': v[1], 'complexo': v[3], 'alias': v[4]}
                              for k, v in FIDUCIAL_DEFS.items()},
            'valve_points': alias_points,
            'valve_times': alias_times,
            'cycle_status': self.cycle_status,
            'corrected_peaks': self.corrected_peaks,
            'ensemble_mean': self.ensemble,
            'ensemble_std': self.ensemble_std,
            'time_axis': self.time_axis,
            'original_peaks': self.peaks,
            'n_cycles': len(self.cycles)
        }
        with open(out_path, 'wb') as f:
            pickle.dump(data, f)
        print(f"   Arquivo: {out_path}")
    
    def run(self):
        plt.show()


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    import sys
    
    print("=" * 60)
    print("FERRAMENTA INTERATIVA DE ANÁLISE SCG")
    print("=" * 60)
    
    # Selecionar paciente
    if len(sys.argv) > 1:
        num = int(sys.argv[1])
    else:
        num = int(input("Número do paciente (1-60): "))
    
    # Tentar carregar do cache primeiro
    CACHE_DIR = DATA_DIR / "_cache_pacientes"
    cache_file = CACHE_DIR / f"paciente_{num}.npz"
    
    if cache_file.exists():
        print(f"\n[Carregando do cache: paciente_{num}.npz]")
        data = np.load(cache_file)
        sinal_clean = data['sinal']
        t = data['t']
        cycles = [data['cycles'][i] for i in range(len(data['cycles']))]
        valid_peaks = [data['peaks'][i] for i in range(len(data['peaks']))]
        fc_est = float(data['fc_est'])
        grupo = str(data['grupo'])
        idade = int(data['idade'])
    else:
        print(f"\n[Cache não encontrado, processando...]")
        # Carregar planilha de metadados (idade/grupo)
        import pandas as pd
        df = pd.read_excel(DATA_DIR / "Pesquisa_PIBIC_sismo.xlsx")
        df.columns = ['Numero', 'Idade', 'Sexo', 'Grupo']
        df['Numero'] = df['Numero'].astype(int)
        df['Grupo'] = df['Grupo'].str.strip().str.lower()
        df['Grupo'] = df['Grupo'].map({'saudável': 'Saudável', 'saudavel': 'Saudável',
                                         'não saudável': 'Não Saudável'})
        df = df.drop_duplicates(subset=['Numero'])

        # Carregar dados brutos
        path = SCG_DIR / f"{num}.xls"
        if not path.exists():
            path = SCG_DIR / f"{num} .xls"
        
        wb = xlrd.open_workbook(str(path))
        ws = wb.sheet_by_name('Raw Data')
        n = ws.nrows - 1
        t = np.array([ws.cell_value(r+1, 0) for r in range(n)])
        acc_mag = np.array([ws.cell_value(r+1, 4) for r in range(n)])
        
        info = df[df['Numero'] == num]
        idade = info['Idade'].values[0]
        grupo = info['Grupo'].values[0]
        
        # Processar
        sinal_clean = remove_noise_sqi(acc_mag, FS)
        fc_est = estimate_fc_autocorr(sinal_clean, FS)
        peaks, fc = detect_initial_beats(sinal_clean, FS)
        cycles, valid_peaks = extract_cycles(sinal_clean, peaks, FS)
        
        # Filtrar ciclos por qualidade
        energies = [np.sum(c**2) for c in cycles]
        med = np.median(energies)
        mad = np.median(np.abs(np.array(energies) - med))
        good = [i for i in range(len(cycles)) if abs(energies[i]-med) < 2*mad]
        cycles = [cycles[i] for i in good]
        valid_peaks = [valid_peaks[i] for i in good]
        
        if len(cycles) < 3:
            print("\n⚠️  Poucos ciclos detectados.")
            sys.exit(1)
    
    print(f"  Paciente #{num} | {grupo} | {idade} anos")
    print(f"  FC estimada: {fc_est:.1f} bpm | Ciclos: {len(cycles)}")
    
    print(f"\n[3] Abrindo ferramenta interativa...")
    tool = SCGInteractiveTool(sinal_clean, t, cycles, valid_peaks, num, grupo, idade)
    tool.run()
    
    print("\n" + "=" * 60)
    print("FERRAMENTA CONCLUÍDA!")
    print("=" * 60)
