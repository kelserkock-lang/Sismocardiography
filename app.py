"""
APP STREAMLIT - ANÁLISE SISMOCARDIOGRÁFICA (SCG)
=================================================
Upload de um arquivo .xls bruto (aba 'Raw Data') e geração automática de:
  - Frequência cardíaca (FC) e variabilidade (HRV: SDNN, RMSSD)
  - Frequência respiratória (FR) estimada
  - Onda média (morfologia) com pontos fiduciais
  - Pontos fiduciais (8 pontos de Sørensen et al., 2018) com tempos e intervalos

Executar:  streamlit run app.py
"""

import io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal
import xlrd
import streamlit as st

FS = 100.0

# ============================================================
# DEFINIÇÕES DOS PONTOS FIDUCIAIS (Sørensen et al., 2018)
# ============================================================
# nome -> (evento, r, cor, complexo, alias válvula)
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
ALIAS_PARA_FIDUCIAL = {v[4]: k for k, v in FIDUCIAL_DEFS.items() if v[4]}
PONTOS = list(FIDUCIAL_DEFS.keys())
INTERVALOS = [
    ('CIV',  'Contração isovolumétrica (MC→AO)', 'Es', 'Gs'),
    ('LVET', 'Ejeção ventricular esq. (AO→AC)',  'Gs', 'Bd'),
    ('RIV',  'Relaxamento isovolumétrico (AC→MO)', 'Bd', 'Fd'),
    ('ENCH', 'Enchimento diastólico (MO→pico E)', 'Fd', 'Gd'),
]

# ============================================================
# VALORES DE REFERÊNCIA (adulto em repouso)
# ============================================================
# métrica -> (mínimo, máximo, unidade, texto de referência, fonte)
VALORES_REF = {
    'fc':    (60, 100, 'bpm', '60–100 bpm',  'Nunan et al., 2010'),
    'fr':    (12, 20,  'rpm', '12–20 rpm',   'Cretikos et al., 2008'),
    'sdnn':  (34, 66,  'ms',  '≈ 50 ms (50 ± 16)',  'Nunan et al., 2010'),
    'rmssd': (27, 57,  'ms',  '≈ 42 ms (42 ± 15)',  'Nunan et al., 2010'),
}
# intervalo -> (mínimo, máximo, texto de referência, fonte)
REF_INTERVALOS = {
    'CIV':  (30, 60,   '30–60 ms',   'Tavakolian, 2016'),
    'LVET': (250, 350, '250–350 ms', 'Sørensen et al., 2018'),
    'RIV':  (60, 100,  '60–100 ms',  'Tavakolian, 2016'),
    'ENCH': (80, 200,  '80–200 ms',  'Sørensen et al., 2018'),
}

# ============================================================
# LITERATURA UTILIZADA
# ============================================================
LITERATURA = [
    'Sørensen K, Schmidt SE, Jensen AS, Søgaard P, Struijk JJ. Definition of '
    'Fiducial Points in the Normal Seismocardiogram. Scientific Reports. '
    '2018;8:15455. doi:10.1038/s41598-018-33675-6.',
    'Nunan D, Sandercock GRH, Brodie DA. A Quantitative Systematic Review of '
    'Normal Values for Short-Term Heart Rate Variability in Healthy Adults. '
    'Pacing and Clinical Electrophysiology. 2010;33(11):1407–1417.',
    'Task Force of the European Society of Cardiology and the North American '
    'Society of Pacing and Electrophysiology. Heart Rate Variability: Standards '
    'of Measurement, Physiological Interpretation, and Clinical Use. '
    'Circulation. 1996;93(5):1043–1065.',
    'Tavakolian K. Systolic Time Intervals and New Measurement Methods. '
    'Cardiovascular Engineering and Technology. 2016;7(2):118–125.',
    'Cretikos MA, Bellomo R, Hillman K, et al. Respiratory rate: the neglected '
    'vital sign. Medical Journal of Australia. 2008;188(11):657–659.',
    'Inan OT, Migeotte PF, Park KS, et al. Ballistocardiography and '
    'Seismocardiography: A Review of Recent Advances. IEEE Journal of '
    'Biomedical and Health Informatics. 2015;19(4):1414–1427.',
]


# ============================================================
# PROCESSAMENTO
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
    b, a = butter_bandpass(2, 39, fs, 4)
    sinal_filt = apply_filter(sinal, b, a)
    threshold = 4 * np.std(sinal_filt)
    mask_bad = np.abs(sinal_filt) > threshold
    sinal_clean = sinal_filt.copy()
    sinal_clean[mask_bad] = np.nan
    if np.any(np.isnan(sinal_clean)):
        x = np.arange(len(sinal_clean))
        good = ~np.isnan(sinal_clean)
        if np.sum(good) > 2:
            sinal_clean = np.interp(x, x[good], sinal_clean[good])
        else:
            sinal_clean = np.nan_to_num(sinal_clean)
    return sinal_clean


def estimate_fc_autocorr(sinal, fs):
    analytic = signal.hilbert(sinal)
    env = np.abs(analytic)
    nyq = fs / 2
    b, a = signal.butter(4, 5 / nyq, btype='low')
    env_filt = signal.filtfilt(b, a, env)
    n = len(env_filt)
    acorr = signal.correlate(env_filt - np.mean(env_filt),
                             env_filt - np.mean(env_filt), mode='full', method='auto')
    acorr = acorr / acorr[n - 1]
    half = acorr[n - 1:]
    half_lags = np.arange(len(half)) / fs
    peaks, _ = signal.find_peaks(half, height=0.05, distance=int(0.15 * fs), prominence=0.03)
    valid = peaks[(half_lags[peaks] >= 0.35) & (half_lags[peaks] <= 1.5)]
    if len(valid) == 0:
        return np.nan
    best = valid[np.argmax(half[valid])]
    return 60.0 / half_lags[best]


def detect_initial_beats(sinal_clean, fs):
    fc = estimate_fc_autocorr(sinal_clean, fs)
    if np.isnan(fc):
        fc = 75.0
    min_dist = max(0.35, 60.0 / fc * 0.7)
    threshold = np.std(sinal_clean) * 0.5
    peaks, _ = signal.find_peaks(sinal_clean, height=threshold,
                                 distance=int(min_dist * fs), prominence=threshold * 0.3)
    return peaks, fc


def extract_cycles(sinal, peaks, fs, before=0.3, after=0.6):
    n_before = int(before * fs)
    n_after = int(after * fs)
    cycles, valid_peaks = [], []
    for p in peaks:
        if p - n_before >= 0 and p + n_after < len(sinal):
            cycles.append(sinal[p - n_before:p + n_after])
            valid_peaks.append(p)
    return cycles, valid_peaks


def calc_hrv(peaks, fs=FS):
    if len(peaks) < 4:
        return np.nan, np.nan, np.nan
    rr = np.diff(np.sort(peaks)) / fs * 1000.0
    rr = rr[(rr > 300) & (rr < 2000)]
    if len(rr) < 3:
        return np.nan, np.nan, np.nan
    fc = 60000.0 / np.mean(rr)
    sdnn = np.std(rr, ddof=1)
    rmssd = np.sqrt(np.mean(np.diff(rr) ** 2))
    return fc, sdnn, rmssd


def estimar_fr(sinal, fs=FS):
    if len(sinal) < int(20 * fs):
        return np.nan
    env = np.abs(signal.hilbert(sinal))
    env = env - np.mean(env)
    b, a = signal.butter(4, 1.0 / (fs / 2), btype='low')
    env = signal.filtfilt(b, a, env)
    f, pxx = signal.welch(env, fs=fs, nperseg=min(len(env), int(30 * fs)))
    banda = (f >= 0.1) & (f <= 0.5)
    if not np.any(banda):
        return np.nan
    return f[banda][np.argmax(pxx[banda])] * 60.0


def detectar_fiduciais(ens, ta):
    s = np.std(ens)
    prom = max(s * 0.10, 1e-6)
    pk, _ = signal.find_peaks(ens, prominence=prom)
    vl, _ = signal.find_peaks(-ens, prominence=prom)
    R = {}
    if len(ens) == 0:
        return R
    win_sis = np.where((ta >= -0.06) & (ta <= 0.12))[0]
    if len(win_sis) == 0:
        win_sis = np.arange(len(ens))
    gs = win_sis[int(np.argmax(ens[win_sis]))]
    R['Gs'] = gs
    esq = pk[pk < gs]
    es = None
    if len(esq) > 0 and (ta[gs] - ta[esq[-1]]) <= 0.16:
        es = int(esq[-1]); R['Es'] = es
    ref = es if es is not None else gs
    esq2 = pk[pk < ref]
    cs = None
    if len(esq2) > 0:
        cs = int(esq2[-1]); R['Cs'] = cs
    ref2 = cs if cs is not None else ref
    esq3 = pk[pk < ref2]
    if len(esq3) > 0:
        R['Bs'] = int(esq3[-1])
    else:
        vesq = vl[vl < ref2]
        if len(vesq) > 0:
            R['Bs'] = int(vesq[-1])
    win_dia = (ta >= 0.13) & (ta <= 0.45)
    dpk = pk[win_dia[pk]]
    bd = None
    if len(dpk) > 0:
        bd = int(dpk[int(np.argmax(ens[dpk]))]); R['Bd'] = bd
    if bd is not None:
        mid = pk[(pk > gs) & (pk < bd)]
    else:
        mid = pk[(pk > gs) & (ta[pk] < 0.15)]
    if len(mid) > 0:
        R['Ks'] = int(mid[int(np.argmax(ens[mid]))])
    fd = None
    if bd is not None:
        dep = pk[pk > bd]
        if len(dep) > 0:
            fd = int(dep[0]); R['Fd'] = fd
    if fd is not None:
        dep2 = pk[pk > fd]
        if len(dep2) > 0:
            R['Gd'] = int(dep2[0])
    return R


def ler_xls(conteudo_bytes):
    """Lê a aba 'Raw Data' de um .xls (colunas: 0=tempo, 4=AccMag)."""
    wb = xlrd.open_workbook(file_contents=conteudo_bytes)
    try:
        ws = wb.sheet_by_name('Raw Data')
    except xlrd.XLRDError:
        ws = wb.sheet_by_index(0)
    n = ws.nrows - 1
    t = np.array([ws.cell_value(r + 1, 0) for r in range(n)], dtype=float)
    acc = np.array([ws.cell_value(r + 1, 4) for r in range(n)], dtype=float)
    return t, acc


def processar(acc_mag):
    """Pipeline completo. Retorna dicionário de resultados."""
    sinal = remove_noise_sqi(acc_mag, FS)
    fc_est = estimate_fc_autocorr(sinal, FS)
    peaks, _ = detect_initial_beats(sinal, FS)
    cycles, valid_peaks = extract_cycles(sinal, peaks, FS)

    # filtrar ciclos por energia
    if len(cycles) > 3:
        energies = np.array([np.sum(c ** 2) for c in cycles])
        med = np.median(energies)
        mad = np.median(np.abs(energies - med)) + 1e-9
        good = [i for i in range(len(cycles)) if abs(energies[i] - med) < 2 * mad]
        cycles = [cycles[i] for i in good]
        valid_peaks = [valid_peaks[i] for i in good]

    fc, sdnn, rmssd = calc_hrv(np.array(peaks))
    if np.isnan(fc):
        fc = fc_est
    fr = estimar_fr(sinal)

    if len(cycles) >= 2:
        ensemble = np.mean(cycles, axis=0)
        ens_std = np.std(cycles, axis=0)
    else:
        ensemble = np.zeros(90)
        ens_std = np.zeros(90)
    ta = np.linspace(-0.3, 0.6, len(ensemble))
    fiduciais = detectar_fiduciais(ensemble, ta)

    return {
        'sinal': sinal, 'peaks': np.array(valid_peaks), 'cycles': cycles,
        'ensemble': ensemble, 'ens_std': ens_std, 'time_axis': ta,
        'fc': fc, 'sdnn': sdnn, 'rmssd': rmssd, 'fr': fr,
        'fiduciais': fiduciais, 'n_ciclos': len(cycles),
    }


# ============================================================
# FIGURAS
# ============================================================

def fig_morfologia(res):
    ens, std, ta = res['ensemble'], res['ens_std'], res['time_axis']
    fid = res['fiduciais']
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(ta, ens, 'b-', lw=2.2, label='Onda média')
    ax.fill_between(ta, ens - std, ens + std, color='blue', alpha=0.12)
    ax.axvline(0, color='red', ls='--', alpha=0.5)
    for nome, (evento, r, cor, cx, alias) in FIDUCIAL_DEFS.items():
        if nome not in fid:
            continue
        idx = fid[nome]
        x, y = ta[idx], ens[idx]
        ax.scatter([x], [y], c=cor, s=130, marker='o',
                   edgecolors='black', linewidths=0.7, zorder=5)
        rot = f'{nome} ({alias})' if alias else nome
        ax.annotate(rot, xy=(x, y), xytext=(0, 11), textcoords='offset points',
                    ha='center', fontsize=9, fontweight='bold', color=cor)
    ax.set_title('Onda Média SCG com Pontos Fiduciais', fontweight='bold')
    ax.set_xlabel('Tempo (s)'); ax.set_ylabel('Aceleração (m/s²)')
    ax.grid(alpha=0.3); ax.legend(loc='upper right', fontsize=8, ncol=2)
    fig.tight_layout()
    return fig


def fig_sinal_fc_fr(res):
    sinal = res['sinal']; peaks = res['peaks']
    t = np.arange(len(sinal)) / FS
    env = np.abs(signal.hilbert(sinal)) - np.mean(np.abs(signal.hilbert(sinal)))
    br, ar = signal.butter(2, [0.1 / (FS / 2), 0.5 / (FS / 2)], btype='band')
    resp = signal.filtfilt(br, ar, env)
    w0, w1 = 10.0, min(30.0, t[-1])
    m = (t >= w0) & (t <= w1)
    bmask = peaks[(t[peaks] >= w0) & (t[peaks] <= w1)] if len(peaks) else np.array([], dtype=int)
    fig, ax1 = plt.subplots(figsize=(10, 4.5))
    ax1.plot(t[m], sinal[m], color='#2c7fb8', lw=0.9, label='SCG filtrado (2–39 Hz)')
    if len(bmask):
        ax1.plot(t[bmask], sinal[bmask], 'v', color='#d7301f', ms=7,
                 label=f'Batimentos → FC = {res["fc"]:.0f} bpm')
    ax1.set_ylabel('Aceleração (m/s²)', color='#2c7fb8')
    ax1.tick_params(axis='y', labelcolor='#2c7fb8')
    ax2 = ax1.twinx()
    denom = np.max(np.abs(resp[m])) + 1e-9
    ax2.plot(t[m], resp[m] / denom, color='#238b45', lw=2.3,
             label=f'Respiração (0,1–0,5 Hz) → FR = {res["fr"]:.0f} rpm')
    ax2.set_ylabel('Respiração (norm.)', color='#238b45')
    ax2.tick_params(axis='y', labelcolor='#238b45'); ax2.set_ylim(-1.6, 1.6)
    l1, la1 = ax1.get_legend_handles_labels()
    l2, la2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, la1 + la2, loc='upper right', fontsize=8)
    ax1.set_title('Captação Simultânea de FC e FR', fontweight='bold')
    ax1.set_xlabel('Tempo (s)'); ax1.grid(alpha=0.25)
    fig.tight_layout()
    return fig


def tabela_fiduciais(res):
    ta = res['time_axis']; fid = res['fiduciais']
    linhas = []
    for nome, (evento, r, cor, cx, alias) in FIDUCIAL_DEFS.items():
        if nome in fid:
            linhas.append({
                'Ponto': nome, 'Alias': alias or '—', 'Evento': evento,
                'Tempo (ms)': f'{ta[fid[nome]]*1000:+.0f}',
                'Complexo': 'Sistólico' if cx == 's' else 'Diastólico',
                'r': r,
            })
    return linhas


def tabela_intervalos(res):
    ta = res['time_axis']; fid = res['fiduciais']
    linhas = []
    for chave, lab, a, b in INTERVALOS:
        ref = REF_INTERVALOS.get(chave, (None, None, '—', ''))[2]
        if a in fid and b in fid:
            dt = (ta[fid[b]] - ta[fid[a]]) * 1000.0
            linhas.append({'Intervalo': chave, 'Descrição': lab,
                           'Valor (ms)': f'{dt:.0f}', 'Referência': ref})
        else:
            linhas.append({'Intervalo': chave, 'Descrição': lab,
                           'Valor (ms)': '—', 'Referência': ref})
    return linhas


def _status(valor, minimo, maximo):
    if valor is None or np.isnan(valor):
        return '—'
    if valor < minimo:
        return 'Abaixo'
    if valor > maximo:
        return 'Acima'
    return 'Normal'


def tabela_referencia(res):
    """Métricas do paciente x faixa de referência da literatura."""
    labels = {'fc': 'FC', 'fr': 'FR', 'sdnn': 'SDNN', 'rmssd': 'RMSSD'}
    linhas = []
    for chave, (mn, mx, un, txt, fonte) in VALORES_REF.items():
        val = res.get(chave, np.nan)
        vtxt = f'{val:.0f} {un}' if val is not None and not np.isnan(val) else '—'
        linhas.append({
            'Métrica': labels[chave], 'Valor': vtxt,
            'Referência': txt, 'Situação': _status(val, mn, mx),
            'Fonte': fonte,
        })
    return linhas



def csv_resultados(res):
    ta = res['time_axis']; fid = res['fiduciais']
    buf = io.StringIO()
    buf.write('metrica,valor\n')
    buf.write(f'FC_bpm,{res["fc"]:.1f}\n')
    buf.write(f'SDNN_ms,{res["sdnn"]:.1f}\n')
    buf.write(f'RMSSD_ms,{res["rmssd"]:.1f}\n')
    buf.write(f'FR_rpm,{res["fr"]:.1f}\n')
    buf.write(f'n_ciclos,{res["n_ciclos"]}\n')
    for nome in PONTOS:
        if nome in fid:
            buf.write(f'tempo_{nome}_ms,{ta[fid[nome]]*1000:.1f}\n')
    for chave, _lab, a, b in INTERVALOS:
        if a in fid and b in fid:
            buf.write(f'{chave}_ms,{(ta[fid[b]]-ta[fid[a]])*1000:.1f}\n')
    return buf.getvalue().encode('utf-8')


def _pagina_tabela(pdf, titulo, colunas, dados, larguras=None, rodape=None):
    from matplotlib.backends.backend_pdf import PdfPages  # noqa: F401
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 retrato
    fig.suptitle(titulo, fontsize=14, fontweight='bold', y=0.96)
    ax = fig.add_axes([0.06, 0.10, 0.88, 0.80]); ax.axis('off')
    tbl = ax.table(cellText=dados, colLabels=colunas, loc='upper center',
                   cellLoc='left', colWidths=larguras)
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.5)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor('#2c3e50'); cell.set_text_props(color='white', fontweight='bold')
        elif r % 2 == 0:
            cell.set_facecolor('#f2f4f6')
    if rodape:
        fig.text(0.06, 0.05, rodape, fontsize=7, color='#555555', wrap=True, va='top')
    pdf.savefig(fig); plt.close(fig)


def gerar_pdf(res, nome_arquivo=''):
    """Gera um relatório PDF completo em memória e devolve bytes."""
    from matplotlib.backends.backend_pdf import PdfPages
    import datetime
    ta = res['time_axis']; fid = res['fiduciais']
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        # --- Capa + resumo ---
        fig = plt.figure(figsize=(8.27, 11.69)); fig.patch.set_facecolor('white')
        fig.text(0.5, 0.90, 'Relatório de Análise Sismocardiográfica',
                 ha='center', fontsize=18, fontweight='bold')
        sub = f'Arquivo: {nome_arquivo}' if nome_arquivo else ''
        data = datetime.datetime.now().strftime('%d/%m/%Y %H:%M')
        fig.text(0.5, 0.865, f'{sub}    Gerado em {data}', ha='center',
                 fontsize=9, color='#555555')
        # cartões de métricas
        cards = [
            ('FC', f'{res["fc"]:.0f} bpm', VALORES_REF['fc'][3]),
            ('FR', f'{res["fr"]:.0f} rpm' if not np.isnan(res['fr']) else '—', VALORES_REF['fr'][3]),
            ('SDNN', f'{res["sdnn"]:.0f} ms' if not np.isnan(res['sdnn']) else '—', VALORES_REF['sdnn'][3]),
            ('RMSSD', f'{res["rmssd"]:.0f} ms' if not np.isnan(res['rmssd']) else '—', VALORES_REF['rmssd'][3]),
            ('Ciclos', f'{res["n_ciclos"]}', 'válidos'),
        ]
        x0, y0, w, h = 0.08, 0.72, 0.40, 0.09
        for i, (lab, val, ref) in enumerate(cards):
            cx = x0 + (i % 2) * 0.44
            cy = y0 - (i // 2) * 0.11
            fig.patches.append(plt.Rectangle((cx, cy), w, h, transform=fig.transFigure,
                                             facecolor='#eef3f7', edgecolor='#2c3e50', lw=1.2))
            fig.text(cx + 0.02, cy + h - 0.025, lab, fontsize=10, color='#2c3e50', fontweight='bold')
            fig.text(cx + 0.02, cy + 0.028, val, fontsize=16, fontweight='bold')
            fig.text(cx + w - 0.02, cy + 0.015, f'ref: {ref}', fontsize=7.5,
                     color='#777777', ha='right')
        fig.text(0.08, 0.16,
                 'Pontos fiduciais segundo Sørensen et al. (2018). Valores de referência '
                 'para adultos em repouso — ver página de referências. Este relatório é '
                 'de apoio à pesquisa e não substitui avaliação clínica.',
                 fontsize=8.5, color='#444444', wrap=True)
        pdf.savefig(fig); plt.close(fig)

        # --- Figuras ---
        f1 = fig_morfologia(res); pdf.savefig(f1); plt.close(f1)
        f2 = fig_sinal_fc_fr(res); pdf.savefig(f2); plt.close(f2)

        # --- Tabela de referência das métricas ---
        ref_rows = [[d['Métrica'], d['Valor'], d['Referência'], d['Situação'], d['Fonte']]
                    for d in tabela_referencia(res)]
        _pagina_tabela(pdf, 'Métricas x Valores de Referência',
                       ['Métrica', 'Valor', 'Referência', 'Situação', 'Fonte'],
                       ref_rows, larguras=[0.14, 0.16, 0.24, 0.16, 0.30])

        # --- Tabela de pontos fiduciais ---
        fid_rows = [[d['Ponto'], d['Alias'], d['Evento'], d['Tempo (ms)'], d['Complexo']]
                    for d in tabela_fiduciais(res)]
        _pagina_tabela(pdf, 'Pontos Fiduciais',
                       ['Ponto', 'Alias', 'Evento', 'Tempo (ms)', 'Complexo'],
                       fid_rows, larguras=[0.10, 0.10, 0.42, 0.16, 0.22])

        # --- Tabela de intervalos ---
        int_rows = [[d['Intervalo'], d['Descrição'], d['Valor (ms)'], d['Referência']]
                    for d in tabela_intervalos(res)]
        _pagina_tabela(pdf, 'Intervalos Temporais',
                       ['Intervalo', 'Descrição', 'Valor (ms)', 'Referência'],
                       int_rows, larguras=[0.14, 0.44, 0.18, 0.24])

        # --- Referências ---
        ref_txt = [[f'{i+1}. {r}'] for i, r in enumerate(LITERATURA)]
        _pagina_tabela(pdf, 'Literatura Utilizada', [''], ref_txt, larguras=[0.98])

    buf.seek(0)
    return buf.getvalue()


# ============================================================
# INTERFACE STREAMLIT
# ============================================================

def run_app():
    st.set_page_config(page_title='Análise SCG', page_icon='🫀', layout='wide')

    st.title('🫀 Análise Sismocardiográfica (SCG)')
    st.caption('Envie um arquivo `.xls` (aba *Raw Data*) para gerar FC, HRV, FR, '
               'morfologia média e pontos fiduciais automaticamente.')

    with st.sidebar:
        st.header('Entrada')
        arquivo = st.file_uploader('Arquivo .xls do paciente', type=['xls'])
        st.markdown('---')
        st.markdown(
            '**Métricas geradas:**\n\n'
            '- FC (bpm) e HRV (SDNN, RMSSD)\n'
            '- FR estimada (rpm)\n'
            '- Onda média + 8 pontos fiduciais\n'
            '- Tempos e intervalos (LVET, CIV, RIV, enchimento)')
        st.caption('Pontos: Sørensen et al., 2018 (Sci. Rep. 8:15455).')

    if arquivo is None:
        st.info('⬅️ Envie um arquivo `.xls` na barra lateral para iniciar a análise.')
        st.stop()

    try:
        with st.spinner('Processando sinal...'):
            _, acc = ler_xls(arquivo.read())
            res = processar(acc)
    except Exception as e:
        st.error(f'Erro ao processar o arquivo: {e}')
        st.stop()

    if res['n_ciclos'] < 2:
        st.warning('Não foi possível extrair ciclos suficientes deste sinal.')
        st.stop()

    # --- Cartões de métricas ---
    st.subheader('Resultados')
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric('FC', f'{res["fc"]:.0f} bpm', help=f'Referência: {VALORES_REF["fc"][3]}')
    c2.metric('FR', f'{res["fr"]:.0f} rpm' if not np.isnan(res['fr']) else '—',
              help=f'Referência: {VALORES_REF["fr"][3]}')
    c3.metric('SDNN', f'{res["sdnn"]:.0f} ms' if not np.isnan(res['sdnn']) else '—',
              help=f'Referência: {VALORES_REF["sdnn"][3]}')
    c4.metric('RMSSD', f'{res["rmssd"]:.0f} ms' if not np.isnan(res['rmssd']) else '—',
              help=f'Referência: {VALORES_REF["rmssd"][3]}')
    c5.metric('Ciclos', f'{res["n_ciclos"]}')

    st.markdown('---')

    # --- Gráficos ---
    g1, g2 = st.columns(2)
    with g1:
        st.markdown('#### Morfologia média + pontos fiduciais')
        st.pyplot(fig_morfologia(res))
    with g2:
        st.markdown('#### Captação de FC e FR')
        st.pyplot(fig_sinal_fc_fr(res))

    st.markdown('---')

    # --- Valores de referência ---
    st.markdown('#### Métricas x valores de referência')
    st.dataframe(tabela_referencia(res), use_container_width=True, hide_index=True)

    st.markdown('---')

    # --- Tabelas ---
    t1, t2 = st.columns(2)
    with t1:
        st.markdown('#### Pontos fiduciais')
        st.dataframe(tabela_fiduciais(res), use_container_width=True, hide_index=True)
    with t2:
        st.markdown('#### Intervalos temporais')
        st.dataframe(tabela_intervalos(res), use_container_width=True, hide_index=True)

    st.markdown('---')

    # --- Literatura ---
    with st.expander('📚 Literatura utilizada'):
        for i, r in enumerate(LITERATURA, 1):
            st.markdown(f'{i}. {r}')

    st.markdown('---')

    # --- Downloads ---
    d1, d2 = st.columns(2)
    with d1:
        st.download_button('⬇️ Baixar resultados (CSV)', data=csv_resultados(res),
                           file_name='resultados_scg.csv', mime='text/csv',
                           use_container_width=True)
    with d2:
        nome = getattr(arquivo, 'name', '')
        pdf_bytes = gerar_pdf(res, nome_arquivo=nome)
        st.download_button('📄 Gerar relatório PDF', data=pdf_bytes,
                           file_name='relatorio_scg.pdf', mime='application/pdf',
                           use_container_width=True)


if __name__ == '__main__':
    run_app()
