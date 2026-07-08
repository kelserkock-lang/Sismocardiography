# 🫀 Análise Sismocardiográfica (SCG) — PIBIC

Ferramenta de pesquisa para **detecção automática de pontos fiduciais** e extração de
métricas cardiorrespiratórias a partir de sinais de **sismocardiografia (SCG)**.

O projeto inclui um **aplicativo web (Streamlit)** que recebe o sinal bruto de um
paciente, executa todo o processamento e apresenta frequência cardíaca, variabilidade
(HRV), frequência respiratória, a morfologia média do ciclo e os 8 pontos fiduciais,
com valores de referência e geração de **relatório em PDF**.

---

## 📱 Como os dados são coletados

Os sinais são adquiridos com um **smartphone** utilizando o aplicativo
**[Phyphox](https://phyphox.org/)**, que dá acesso aos sensores do aparelho:

- O smartphone é posicionado sobre o **tórax** do participante (região do esterno/precórdio),
  com a pessoa em **repouso**, deitada em decúbito dorsal.
- É utilizado o **acelerômetro** do celular, que capta as micro-vibrações da parede torácica
  geradas pela atividade mecânica do coração (abertura/fechamento de válvulas, ejeção e
  enchimento ventricular) — o princípio da sismocardiografia.
- Os dados são exportados pelo Phyphox em planilha **`.xls`** (aba *Raw Data*), contendo o
  tempo e as componentes de aceleração. A análise utiliza a **magnitude da aceleração**
  ($\sqrt{a_x^2 + a_y^2 + a_z^2}$), amostrada a **100 Hz**.

> ⚠️ **Aviso:** ferramenta de apoio à pesquisa (PIBIC). Não é um dispositivo médico e não
> substitui avaliação clínica ou exames como ECG/ecocardiograma.

---

## 🔬 O que a ferramenta calcula

| Métrica | Descrição |
|---|---|
| **FC** | Frequência cardíaca (bpm), via autocorrelação do envelope + detecção de batimentos |
| **HRV** | Variabilidade da FC: **SDNN** e **RMSSD** (ms) |
| **FR** | Frequência respiratória (rpm), pelo envelope respiratório (Hilbert + banda 0,1–0,5 Hz) |
| **Morfologia** | Onda média do ciclo cardíaco (ensemble average) ± desvio |
| **Pontos fiduciais** | 8 pontos (Sørensen et al., 2018): MC, AO, AC, MO e picos associados |
| **Intervalos** | CIV, LVET (ejeção), RIV e enchimento diastólico (ms) |

Os pontos fiduciais seguem a nomenclatura de **Sørensen et al. (2018)**, com destaque para
os eventos valvares: **MC** (fechamento mitral), **AO** (abertura aórtica),
**AC** (fechamento aórtico) e **MO** (abertura mitral).

---

## 🗂️ Estrutura do projeto

| Arquivo | Função |
|---|---|
| `app.py` | Aplicativo web Streamlit (upload do `.xls` → análise → PDF) |
| `ferramenta_interativa.py` | Ferramenta desktop para marcação manual dos pontos fiduciais |
| `identificar_pontos_automatico.py` | Detecção automática dos pontos em lote (todos os pacientes) |
| `gerar_relatorio_word.py` | Relatório Word comparando grupos (saudável × não saudável) |
| `GUIA_PONTOS_FIDUCIAIS_SCG.md` | Guia técnico sobre os pontos fiduciais e sua marcação |
| `requirements.txt` | Dependências do projeto |

---

## ▶️ Executando localmente

```bash
# 1. (recomendado) criar um ambiente virtual
python -m venv .venv
# Windows
.venv\Scripts\activate
# Linux/Mac
source .venv/bin/activate

# 2. instalar as dependências
pip install -r requirements.txt

# 3. iniciar o aplicativo
streamlit run app.py
```

O app abre em `http://localhost:8501`. Envie um arquivo `.xls` (exportado do Phyphox) pela
barra lateral para iniciar a análise.

---

## 📚 Literatura utilizada

1. Sørensen K, Schmidt SE, Jensen AS, Søgaard P, Struijk JJ. *Definition of Fiducial Points
   in the Normal Seismocardiogram.* Scientific Reports. 2018;8:15455.
2. Nunan D, Sandercock GRH, Brodie DA. *A Quantitative Systematic Review of Normal Values for
   Short-Term Heart Rate Variability in Healthy Adults.* PACE. 2010;33(11):1407–1417.
3. Task Force of the ESC/NASPE. *Heart Rate Variability: Standards of Measurement,
   Physiological Interpretation, and Clinical Use.* Circulation. 1996;93(5):1043–1065.
4. Tavakolian K. *Systolic Time Intervals and New Measurement Methods.* Cardiovascular
   Engineering and Technology. 2016;7(2):118–125.
5. Cretikos MA, et al. *Respiratory rate: the neglected vital sign.* MJA. 2008;188(11):657–659.
6. Inan OT, et al. *Ballistocardiography and Seismocardiography: A Review of Recent Advances.*
   IEEE J. Biomed. Health Inform. 2015;19(4):1414–1427.

---

## 🔒 Privacidade dos dados

Os dados de pacientes (**arquivos `.xls`, `.npz` e planilhas**) **não são versionados** neste
repositório (ver `.gitignore`). O aplicativo processa o arquivo enviado **em memória** e não o
armazena. Ao utilizar dados de participantes, garanta **anonimização** e **consentimento
informado**, conforme as normas éticas de pesquisa.

---

*Projeto desenvolvido no âmbito do PIBIC — Iniciação Científica.*
