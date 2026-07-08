# Pontos Fiduciais da Sismocardiografia (SCG) — Guia de Referência

> Documento-base consolidando o conhecimento sobre identificação e marcação de
> pontos fiduciais em sismocardiografia, sua relação com eventos fisiológicos,
> e a importância para o template médio, contagem de FC e demais análises.
>
> **Fonte principal:** Sørensen, K., Schmidt, S. E., Jensen, A. S., Søgaard, P.,
> & Struijk, J. J. (2018). *Definition of Fiducial Points in the Normal
> Seismocardiogram.* Scientific Reports, 8:15455.
> DOI: 10.1038/s41598-018-33675-6

---

## 1. O que é sismocardiografia (SCG)

A sismocardiografia mede as **vibrações mecânicas produzidas pelo coração**,
captadas na parede torácica por um **acelerômetro** (tipicamente sobre o
processo xifoide). Diferente do ECG (elétrico), o SCG reflete os **eventos
mecânicos** do ciclo cardíaco: abertura/fechamento de válvulas, ejeção e
enchimento ventricular.

- Sinal **determinístico** e reprodutível entre sujeitos.
- Mais detalhado (mais ondas) que o ECG.
- Tecnologia viável hoje graças a acelerômetros MEMS (pequenos, leves,
  sensíveis, vestíveis).
- Aplicações: timing de eventos cardíacos, PEP/LVET, tempo de trânsito de
  pulso (PTT), monitoramento de insuficiência cardíaca.

---

## 2. Esquema de rotulagem dos pontos fiduciais (A–M)

Os pontos fiduciais são os **picos e vales** característicos da onda SCG média.
Sørensen et al. rotularam-nos com letras **A a M**, separadas em dois complexos:

- **Complexo sistólico** — subscrito `s` (ex.: `Aₛ`, `Fₛ`, `Gₛ`)
- **Complexo diastólico** — subscrito `d` (ex.: `A_d`, `B_d`, `E_d`)

### 2.1 Complexo sistólico (após o pico R do ECG, t = 0)

A rotulagem parte da **primeira deflexão negativa clara** do complexo principal
seguindo o pico R:

| Ponto | Definição morfológica |
|-------|-----------------------|
| `Fₛ`  | Primeira deflexão negativa clara do complexo principal (âncora) |
| `Eₛ`  | Pequeno pico positivo **antes** de F (às vezes só um "ombro") |
| `Dₛ`  | Pico positivo antes de E (pula a deflexão negativa entre E e D) |
| `Cₛ`  | Vale antes de D |
| `Bₛ`  | Pico antes de C |
| `Aₛ`  | Primeiro pico (frequentemente **ausente**) |
| `Gₛ`  | Pico entre F e I (abertura aórtica) |
| `Hₛ`  | Segundo pico entre F e I (só se houver > 15 ms de Gₛ) |
| `Iₛ`  | Vale após o pico mais alto do complexo sistólico |
| `Jₛ`  | Pico após Iₛ |
| `Kₛ`  | Vale/pico antes de Mₛ (pico de ejeção sistólica) |
| `Lₛ`  | Pico antes de Mₛ |
| `Mₛ`  | Último vale marcando o fim da sístole |

**Regras práticas de Gₛ/Hₛ:** entre `Fₛ` e `Iₛ` podem existir 1, 2 ou 3 picos.
- 1 pico → `Gₛ`.
- 2 picos → `Gₛ` e `Hₛ`, **mas só se houver > 15 ms entre eles**; caso
  contrário, apenas `Gₛ`.
- 3 picos → considera-se apenas os 2 mais próximos de `Fₛ`.

### 2.2 Complexo diastólico (alinhado à 2ª bulha, S2)

| Ponto | Definição morfológica |
|-------|-----------------------|
| `A_d` | Primeiro vale (ou início da aceleração positiva diastólica) |
| `B_d` | Pico seguinte (às vezes um "ombro" de aceleração negativa) |
| `C_d` | Vale íngreme após B |
| `D_d` | Pico após C |
| `E_d` | Vale mais profundo após D |
| `F_d`, `G_d`, `H_d`, `I_d` | Picos e vales subsequentes |

---

## 3. Os 8 pontos fiduciais com correlação fisiológica validada

De todos os pontos A–M, **apenas 8** têm correlação estatisticamente
significativa com eventos cardíacos (validados por ultrassom: Doppler pulsado e
Tissue Doppler). **Estes são os pontos que importam na prática clínica:**

| Ponto | Evento fisiológico | Diferença média (ms) | Correlação r | Alias |
|-------|--------------------|----------------------|--------------|-------|
| `Bₛ`  | Sístole atrial (início da onda A) | −2 (±16) | 0,75 | — |
| `Cₛ`  | Pico de influxo atrial | +13 (±19) | 0,63 | — |
| `Eₛ`  | **Fechamento da válvula mitral** | +4 (±11) | 0,71 | **MC** |
| `Gₛ`  | **Abertura da válvula aórtica** | −3 (±11) | 0,60 | **AO** |
| `Kₛ`  | Pico de ejeção sistólica | +13 (±23) | 0,42 | — |
| `B_d` | **Fechamento da válvula aórtica** | −5 (±12) | **0,94** | **AC** |
| `F_d` | **Abertura da válvula mitral** | −7 (±19) | 0,87 | **MO** |
| `G_d` | Pico de enchimento ventricular (onda E) | −18 (±28) | 0,79 | — |

> **Nota sobre OCR:** o artigo original tem erros de digitalização ("Sₐ", "Uₛ")
> que correspondem, respectivamente, aos pontos diastólicos `F_d` (MO) e `G_d`.

### 3.1 Ordem de confiabilidade
1. **`B_d` (AC — fechamento aórtico):** o mais robusto (**r = 0,94**). Pico
   proeminente no complexo diastólico.
2. **`F_d` (MO — abertura mitral):** r = 0,87.
3. **`G_d` (pico da onda E):** r = 0,79.
4. **`Bₛ` (sístole atrial):** r = 0,75.
5. **`Eₛ` (MC), `Cₛ`, `Gₛ` (AO):** r = 0,60–0,71.
6. **`Kₛ` (ejeção):** o mais fraco (r = 0,42).

> **Ponto de partida recomendado na marcação:** `Gₛ` (AO) é o pico sistólico
> dominante — o mais fácil de localizar visualmente. Depois `B_d` (AC), o mais
> confiável. Os demais se posicionam por ordenação relativa a esses dois.

---

## 4. Marcação adequada — metodologia

### 4.1 Pré-processamento do sinal (do artigo)
1. Filtro passa-baixa Butterworth 1ª ordem, corte **90 Hz** (forward-backward).
2. Filtro passa-alta Butterworth 3ª ordem, corte **0,05 Hz** (forward-backward).
3. Segmentação em batimentos pelo **pico R do ECG**.
4. Cálculo do **batimento médio (ensemble)**, com sístole e diástole alinhadas
   e promediadas separadamente (para respeitar a duração variável da sístole).
5. Rejeição de batimentos ruidosos antes de promediar.

### 4.2 Duplo alinhamento (crucial)
O batimento médio é calculado com **dois alinhamentos distintos**:

- **Complexo sistólico:** alinhado ao **pico R do ECG** (t = 0).
- **Complexo diastólico:** alinhado à **2ª bulha cardíaca (S2)**, localizada por
  autocorrelação do envelope do som cardíaco (banda 50–500 Hz do IC4).

**Por quê?** A variabilidade do intervalo S1–S2 faz o alinhamento único (pelo R)
"borrar" os eventos diastólicos. O alinhamento por S2 garante que nenhum evento
diastólico seja perdido na promediação.

### 4.3 Processo de anotação
- Anotação **manual**, iterativa, com o operador percorrendo os batimentos
  médios várias vezes.
- Passo 1: marcar **todos** os picos e vales significativos.
- Passo 2: rotular os eventos (A–M) de forma **consistente entre sujeitos** —
  pontos de mesma característica recebem a mesma letra.
- O operador deve estar **cego** para as imagens de ultrassom (evita viés).
- Nem todos os pontos aparecem em todos os sujeitos (ex.: `Aₛ` presente em só
  ~26% dos sinais; `Hₛ` em ~24%).

### 4.4 Adaptação sem ECG (importante para outros projetos)
Se **não houver ECG** (só o acelerômetro), como no nosso pipeline:
- O ensemble é alinhado ao **pico de aceleração dominante** (não ao pico R).
- Nesse caso, **t = 0 corresponde ao pico sistólico dominante**, não ao R.
- Consequência: os **pontos sistólicos** são marcáveis diretamente e de forma
  confiável; os **pontos diastólicos** (`B_d`, `F_d`, `G_d`) aparecem, mas o
  alinhamento **não é idêntico** ao do artigo (falta a âncora de S2).
- Os tempos absolutos passam a ser relativos ao pico sistólico; por isso,
  **prefira intervalos entre pontos** (ver §6) a tempos absolutos.

---

## 5. Importância da identificação dos pontos fiduciais

### 5.1 Para entender o template (onda média)
- O **template médio (ensemble)** só faz sentido fisiológico quando ancorado em
  pontos fiduciais consistentes. Sem eles, a onda média é apenas uma forma;
  com eles, cada deflexão ganha **significado mecânico** (abertura/fechamento
  valvar, ejeção, enchimento).
- A promediação **reduz o ruído** e revela a morfologia reprodutível — mas exige
  **alinhamento correto** (senão os eventos se cancelam/borram).
- Pontos fiduciais permitem **comparar morfologias entre grupos** (ex.:
  saudável × não saudável) de forma objetiva, ponto a ponto.
- Servem de **âncora para normalização temporal** entre sujeitos com FC
  diferentes.

### 5.2 Para contar corretamente a FC
- A FC é derivada da **detecção de batimentos** (picos) no sinal. Conhecer a
  morfologia (qual é o pico dominante = `Gₛ`/AO) evita erros de contagem.
- **Erro comum: dupla contagem.** O complexo sistólico tem vários picos
  próximos (`Fₛ`–`Gₛ`–`Iₛ`–`Jₛ`). Um detector ingênuo conta 2+ picos por
  batimento → **FC superestimada** (ex.: 76 bpm virando 110 bpm).
  - **Solução:** impor **distância mínima** entre picos (ex.:
    `distance ≈ 0.4 s`, ~150 bpm máx.) e usar **proeminência** para pegar só
    o pico dominante do complexo.
- **Erro oposto: subcontagem.** Distância mínima grande demais funde batimentos
  próximos → FC subestimada.
- **Boa prática:** estimar a FC preliminar por **autocorrelação do envelope**
  (banda ~5 Hz) e usar isso para ajustar a distância mínima da detecção final.
- Validar a FC por **regularidade dos intervalos RR** (CV baixo = detecção
  limpa). Filtrar RR fisiologicamente implausíveis (< 300 ms ou > 2000 ms).

### 5.3 Para variabilidade da FC (HRV)
- A HRV depende de **intervalos RR precisos** entre batimentos → depende da
  detecção correta do **mesmo ponto fiducial** em cada ciclo.
- Métricas no domínio do tempo:
  - **SDNN** = desvio-padrão dos intervalos RR (ms).
  - **RMSSD** = raiz quadrada da média dos quadrados das diferenças sucessivas
    de RR (ms) — reflete atividade parassimpática.
- Marcar sempre o **mesmo ponto** (ex.: `Gₛ`/AO) em todos os ciclos é essencial
  para não introduzir jitter artificial na HRV.

---

## 6. Intervalos temporais entre pontos (com significado clínico)

Mais robustos que tempos absolutos (independem do alinhamento):

| Intervalo | Pontos | Significado |
|-----------|--------|-------------|
| **CIV** (contração isovolumétrica) | `Eₛ` → `Gₛ` (MC → AO) | Tempo entre fechamento mitral e abertura aórtica |
| **LVET** (tempo de ejeção do VE) | `Gₛ` → `B_d` (AO → AC) | Duração da ejeção ventricular esquerda |
| **RIV** (relaxamento isovolumétrico) | `B_d` → `F_d` (AC → MO) | Tempo entre fechamento aórtico e abertura mitral |
| **Enchimento** | `F_d` → `G_d` (MO → pico E) | Início do enchimento ventricular rápido |

- **PEP** (período pré-ejeção) e **LVET** são clássicos indicadores de
  contratilidade e de status hemodinâmico.
- Diferenças de LVET/RIV entre grupos podem indicar alterações da função
  sistólica/diastólica.

---

## 7. Extração da frequência respiratória (FR) do próprio SCG

O SCG carrega **modulação respiratória** (a amplitude dos batimentos varia com a
respiração). É possível extrair a FR de um **único acelerômetro**:

1. Filtrar o SCG na **banda cardíaca 2–39 Hz** (passa-banda Butterworth).
2. Calcular o **envelope de amplitude** via transformada de **Hilbert**.
3. Filtrar o envelope na **banda respiratória 0,1–0,5 Hz** (6–30 rpm).
4. **FR** = frequência do pico dominante nessa banda × 60 (rpm).

- O **espectro do envelope** mostra **dois picos**: um respiratório (~0,2–0,4 Hz)
  e um cardíaco (~1–1,7 Hz) — evidenciando que um único sensor capta FC **e** FR.
- Sem sinal respiratório dedicado, a FR é uma **estimativa** (não padrão-ouro).

---

## 8. Bandas de filtragem — resumo

| Objetivo | Filtro | Faixa |
|----------|--------|-------|
| Componente cardíaca (SCG) | Passa-banda Butterworth | **2–39 Hz** |
| Pré-proc. do artigo | Passa-baixa 90 Hz + passa-alta 0,05 Hz | — |
| Som cardíaco (S1/S2) | Passa-banda | 50–500 Hz |
| Envelope respiratório | Passa-banda sobre o envelope | **0,1–0,5 Hz** |
| Estimativa de FC (envelope) | Passa-baixa | ~5 Hz |

> Sempre usar filtragem **forward-backward** (`filtfilt`) para **fase zero**
> (não deslocar os pontos fiduciais no tempo). Aplicar **padding** (reflexão)
> nas bordas para evitar artefatos de transição.

---

## 9. Boas práticas e armadilhas (aprendizados)

- **Alinhamento correto é tudo:** sistólico pelo R (ou pico dominante),
  diastólico por S2. Alinhamento único borra a diástole.
- **Rejeitar batimentos ruidosos** antes de promediar o ensemble.
- **Distância mínima + proeminência** na detecção de picos evitam dupla
  contagem da FC.
- **Filtragem de fase zero** (`filtfilt`) preserva a posição temporal dos
  pontos — filtros causais deslocariam os fiduciais.
- **Preferir intervalos entre pontos** a tempos absolutos quando não há ECG.
- Nem todo ponto existe em todo sujeito — o esquema é **flexível** (marca-se o
  que está presente).
- **Padronizar o ponto de referência** (ex.: sempre `Gₛ`/AO) entre ciclos e
  entre sujeitos para comparabilidade e HRV consistente.
- Validar a FC pela **regularidade dos RR** (CV baixo indica sinal limpo).
- Para escolher um **sinal representativo/limpo**, usar uma métrica de **SNR**:
  amplitude do batimento médio ÷ dispersão entre ciclos.
- A correlação `B_d` (AC) = 0,94 torna-o a **âncora diastólica mais confiável**;
  `Gₛ` (AO) é a **âncora sistólica** visualmente mais evidente.

---

## 10. Nomenclatura rápida (colar em código)

```text
# Pontos validados: nome -> (evento, r, complexo, alias de válvula)
Bs -> Sístole atrial            r=0.75  sistólico   (—)
Cs -> Pico de influxo atrial    r=0.63  sistólico   (—)
Es -> Fechamento mitral         r=0.71  sistólico   (MC)
Gs -> Abertura aórtica          r=0.60  sistólico   (AO)
Ks -> Pico de ejeção sistólica  r=0.42  sistólico   (—)
Bd -> Fechamento aórtico        r=0.94  diastólico  (AC)
Fd -> Abertura mitral           r=0.87  diastólico  (MO)
Gd -> Pico da onda E            r=0.79  diastólico  (—)

# Intervalos:
CIV  = t(Gs) - t(Es)   # MC -> AO
LVET = t(Bd) - t(Gs)   # AO -> AC
RIV  = t(Fd) - t(Bd)   # AC -> MO
ENCH = t(Gd) - t(Fd)   # MO -> pico E
```

---

## 11. Referência

Sørensen, K., Schmidt, S. E., Jensen, A. S., Søgaard, P., & Struijk, J. J.
(2018). **Definition of Fiducial Points in the Normal Seismocardiogram.**
*Scientific Reports*, 8(1), 15455. https://doi.org/10.1038/s41598-018-33675-6

- 45 sujeitos saudáveis; 42 sinais analisados.
- SCG + ECG a 5000 Hz; ultrassom (Doppler pulsado + Tissue Doppler) como
  referência dos eventos.
- 8 eventos fisiológicos ↔ pontos fiduciais reprodutíveis e bem definidos.
