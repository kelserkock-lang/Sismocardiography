"""
IDENTIFICAÇÃO AUTOMÁTICA DE PONTOS FIDUCIAIS SCG
=================================================
Detecta automaticamente os 8 pontos fiduciais validados por
Sørensen et al. (2018, Sci. Rep. 8:15455) na onda média (ensemble)
de cada paciente do cache, salva as anotações (formato compatível
com a ferramenta interativa) e gera gráficos de conferência (QC).

A onda média está alinhada pelo pico sistólico dominante (t ~ 0 s).
A detecção é heurística, baseada na ordenação de picos/vales em
relação a esse pico e ao complexo diastólico de fechamento aórtico.
Os pontos gerados são ESTIMATIVAS iniciais que podem ser revisadas
manualmente na ferramenta_interativa.py.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # backend não interativo (processamento em lote)
import matplotlib.pyplot as plt
from scipy import signal
from pathlib import Path
import pickle
import csv
import warnings
warnings.filterwarnings('ignore')

DATA_DIR = Path(r"C:\Users\kelse\OneDrive\Área de Trabalho\DADOS_PIBIC")
CACHE_DIR = DATA_DIR / "_cache_pacientes"
OUT_DIR = DATA_DIR / "pontos_automaticos"
FS = 100.0

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
ALIAS_PARA_FIDUCIAL = {v[4]: k for k, v in FIDUCIAL_DEFS.items() if v[4]}


# ============================================================
# DETECÇÃO AUTOMÁTICA
# ============================================================

def detectar_fiduciais(ens, ta, fs=FS):
    """
    Detecta os 8 pontos fiduciais na onda média.

    Parâmetros
    ----------
    ens : np.ndarray  -> onda média (ensemble)
    ta  : np.ndarray  -> eixo de tempo correspondente (s)

    Retorna
    -------
    dict {nome_ponto: indice_no_ensemble}
    """
    s = np.std(ens)
    prom = max(s * 0.10, 1e-6)

    pk, _ = signal.find_peaks(ens, prominence=prom)
    vl, _ = signal.find_peaks(-ens, prominence=prom)
    R = {}

    if len(ens) == 0:
        return R

    # ---- Âncora: pico sistólico dominante (Gs / AO) ----
    win_sis = np.where((ta >= -0.06) & (ta <= 0.12))[0]
    if len(win_sis) == 0:
        win_sis = np.arange(len(ens))
    gs = win_sis[int(np.argmax(ens[win_sis]))]
    R['Gs'] = gs

    # ---- Es (MC): pico positivo imediatamente antes de Gs ----
    esq = pk[pk < gs]
    es = None
    if len(esq) > 0 and (ta[gs] - ta[esq[-1]]) <= 0.16:
        es = int(esq[-1])
        R['Es'] = es

    # ---- Cs (pico de influxo atrial): pico antes de Es (ou de Gs) ----
    ref = es if es is not None else gs
    esq2 = pk[pk < ref]
    cs = None
    if len(esq2) > 0:
        cs = int(esq2[-1])
        R['Cs'] = cs

    # ---- Bs (sístole atrial): pico/onset antes de Cs ----
    ref2 = cs if cs is not None else ref
    esq3 = pk[pk < ref2]
    if len(esq3) > 0:
        R['Bs'] = int(esq3[-1])
    else:
        # fallback: vale (onset) mais próximo à esquerda na região atrial
        vesq = vl[vl < ref2]
        if len(vesq) > 0:
            R['Bs'] = int(vesq[-1])

    # ---- Bd (AC): pico positivo mais proeminente na diástole [0.13, 0.45] ----
    win_dia = (ta >= 0.13) & (ta <= 0.45)
    dpk = pk[win_dia[pk]]
    bd = None
    if len(dpk) > 0:
        bd = int(dpk[int(np.argmax(ens[dpk]))])
        R['Bd'] = bd

    # ---- Ks (pico de ejeção sistólica): pico entre Gs e Bd ----
    if bd is not None:
        mid = pk[(pk > gs) & (pk < bd)]
    else:
        mid = pk[(pk > gs) & (ta[pk] < 0.15)]
    if len(mid) > 0:
        R['Ks'] = int(mid[int(np.argmax(ens[mid]))])

    # ---- Fd (MO): próximo pico proeminente após Bd ----
    fd = None
    if bd is not None:
        dep = pk[pk > bd]
        if len(dep) > 0:
            fd = int(dep[0])
            R['Fd'] = fd

    # ---- Gd (pico onda E): próximo pico proeminente após Fd ----
    if fd is not None:
        dep2 = pk[pk > fd]
        if len(dep2) > 0:
            R['Gd'] = int(dep2[0])

    return R


# ============================================================
# CARGA / PLOT
# ============================================================

def carregar_ensemble(num):
    """Carrega onda média, eixo de tempo e metadados do cache."""
    cache_file = CACHE_DIR / f"paciente_{num}.npz"
    if not cache_file.exists():
        return None
    d = np.load(cache_file)
    cycles = [d['cycles'][i] for i in range(len(d['cycles']))]
    if len(cycles) < 2:
        return None
    ens = np.mean(cycles, axis=0)
    std = np.std(cycles, axis=0)
    ta = np.linspace(-0.3, 0.6, len(ens))
    grupo = str(d['grupo']) if 'grupo' in d else '?'
    idade = int(d['idade']) if 'idade' in d else -1
    peaks = d['peaks'] if 'peaks' in d else np.array([])
    return {
        'ensemble': ens, 'std': std, 'time_axis': ta,
        'grupo': grupo, 'idade': idade, 'n_cycles': len(cycles),
        'peaks': peaks,
    }


def gerar_qc_plot(num, dados, pontos, out_path):
    """Gera figura de conferência com os pontos detectados."""
    ens = dados['ensemble']
    std = dados['std']
    ta = dados['time_axis']

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.plot(ta, ens, 'b-', lw=2.2, label='Onda média')
    ax.fill_between(ta, ens - std, ens + std, color='blue', alpha=0.12)
    ax.axvline(0, color='red', ls='--', alpha=0.5, lw=1.2)
    ax.axhline(0, color='gray', ls=':', alpha=0.4)

    for nome, (evento, r, cor, cx, alias) in FIDUCIAL_DEFS.items():
        if nome not in pontos:
            continue
        idx = pontos[nome]
        x, y = ta[idx], ens[idx]
        ax.scatter([x], [y], c=cor, s=140, marker='o',
                   edgecolors='black', linewidths=0.8, zorder=5)
        rot = f'{nome} ({alias})' if alias else nome
        ax.annotate(rot, xy=(x, y), xytext=(0, 12), textcoords='offset points',
                    ha='center', fontsize=10, fontweight='bold', color=cor, zorder=6)

    ax.set_title(f'Paciente #{num} | {dados["grupo"]} | {dados["idade"]} anos | '
                 f'n={dados["n_cycles"]} ciclos  —  Pontos fiduciais (AUTO)',
                 fontsize=13, fontweight='bold')
    ax.set_xlabel('Tempo (s)')
    ax.set_ylabel('Aceleração (m/s²)')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def salvar_anotacoes(num, dados, pontos):
    """Salva pkl compatível com a ferramenta interativa."""
    ta = dados['time_axis']
    alias_points, alias_times = {}, {}
    for alias, nome in ALIAS_PARA_FIDUCIAL.items():
        if nome in pontos:
            idx = pontos[nome]
            alias_points[alias] = idx
            alias_times[alias] = float(ta[idx])
    data = {
        'patient': num,
        'grupo': dados['grupo'],
        'idade': dados['idade'],
        'metodo': 'automatico',
        'fiducial_points': pontos,
        'fiducial_times': {k: float(ta[v]) for k, v in pontos.items()},
        'fiducial_defs': {k: {'evento': v[0], 'r': v[1], 'complexo': v[3], 'alias': v[4]}
                          for k, v in FIDUCIAL_DEFS.items()},
        'valve_points': alias_points,
        'valve_times': alias_times,
        'ensemble_mean': dados['ensemble'],
        'ensemble_std': dados['std'],
        'time_axis': ta,
        'original_peaks': dados['peaks'],
        'n_cycles': dados['n_cycles'],
    }
    out_path = DATA_DIR / f"anotacoes_auto_paciente_{num}.pkl"
    with open(out_path, 'wb') as f:
        pickle.dump(data, f)
    return out_path


# ============================================================
# MAIN
# ============================================================

def main():
    OUT_DIR.mkdir(exist_ok=True)
    arquivos = sorted(CACHE_DIR.glob("paciente_*.npz"),
                      key=lambda p: int(p.stem.split('_')[1]))
    if not arquivos:
        print("Nenhum paciente encontrado no cache.")
        return

    print("=" * 64)
    print("IDENTIFICAÇÃO AUTOMÁTICA DE PONTOS FIDUCIAIS SCG")
    print("=" * 64)
    print(f"Pacientes no cache: {len(arquivos)}\n")

    resumo = []
    nomes = list(FIDUCIAL_DEFS.keys())
    ok, falhas = 0, 0

    for arq in arquivos:
        num = int(arq.stem.split('_')[1])
        dados = carregar_ensemble(num)
        if dados is None:
            print(f"  #{num:>3}: sem dados suficientes — pulado")
            falhas += 1
            continue

        pontos = detectar_fiduciais(dados['ensemble'], dados['time_axis'])
        salvar_anotacoes(num, dados, pontos)
        gerar_qc_plot(num, dados, pontos, OUT_DIR / f"paciente_{num}.png")

        ta = dados['time_axis']
        linha = {'paciente': num, 'grupo': dados['grupo'],
                 'idade': dados['idade'], 'n_ciclos': dados['n_cycles'],
                 'n_pontos': len(pontos)}
        for nm in nomes:
            linha[f't_{nm}'] = round(float(ta[pontos[nm]]), 3) if nm in pontos else ''
        resumo.append(linha)

        det = " ".join(nomes[i] if nomes[i] in pontos else '·' for i in range(len(nomes)))
        print(f"  #{num:>3}: {len(pontos)}/8 pontos  [{det}]")
        ok += 1

    # Resumo CSV
    csv_path = DATA_DIR / "pontos_automaticos_resumo.csv"
    if resumo:
        campos = ['paciente', 'grupo', 'idade', 'n_ciclos', 'n_pontos'] + [f't_{n}' for n in nomes]
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=campos)
            w.writeheader()
            w.writerows(resumo)

    print("\n" + "=" * 64)
    print(f"CONCLUÍDO: {ok} processados, {falhas} pulados")
    print(f"  Gráficos QC:  {OUT_DIR}")
    print(f"  Anotações:    anotacoes_auto_paciente_N.pkl")
    print(f"  Resumo:       {csv_path.name}")
    print("=" * 64)


if __name__ == '__main__':
    main()
