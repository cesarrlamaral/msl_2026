# Notas editoriais — paper_msl_2026.docx

Pendências e decisões editoriais levantadas durante a reorganização, verificação de referências e refinamento científico (2026-09-05). Nada disto vive no corpo do manuscrito — é o lugar certo pra esse tipo de nota operacional (mesmo padrão do `CITATIONS_TO_VERIFY.md` do hydrovent_field).

## Refinamento científico #2 (2026-09-05) — busca ampla na literatura
Busca web ampla pediu achados relevantes pra pesquisa. Três referências novas reais, verificadas com fonte primária, incorporadas ao paper:
- **Ehresmann et al. (2023, Icarus 393:115035)** — dados de nov/2019-out/2020 (mínimo solar), ~50% de aumento de dose entre 2012-13 e 2019-20. Confirma independentemente (e estende) nosso próprio achado Period1→2.
- **Guo et al. (2021, A&A Review 29:8)** — revisão completa até dez/2020, mesmo ~50%.
- **Zhang et al. (2026, Space Weather)** — limites de missão tripulada; nuance de que trânsito em máximo solar pode ter risco total menor mesmo com dose de superfície maior no mínimo (efeito de SEP em trânsito), usada na Conclusion.

**Erro de busca detectado e corrigido antes de entrar no paper**: o resumo automático da primeira busca atribuiu ao Guo (2015) um número específico ("~25% da dose, 55,5 µGy/dia de diferença sazonal") que **não existe no PDF real** (já estava na nossa bibliografia) — o "25%" real ali é sobre a variação sazonal da PRESSÃO, não da dose. Fui direto no PDF confirmar antes de escrever qualquer coisa, e achei a frase real e muito mais útil: Guo (2015) diz que "the seasonal pressure influence is much smaller than the longer-term effects driven by solar modulation" — isso bate exatamente com nosso próprio achado negativo (wavelet/decomposição não resolve periodicidade sazonal), e agora está citado corretamente no texto (Seção 3.3).

**Reposicionamento da contribuição do paper**: como Guo (2021) e Ehresmann (2023) já estenderam o registro até o mínimo solar de verdade (com modelagem mais completa), a Introduction e a Conclusion foram reescritas pra não reivindicar mais "mais dados que a literatura" como diferencial — a contribuição real e defensável é o **teste espectral direto** (wavelet + decomposição, ambos com escala corrigida) da expectativa sazonal-vs-tendência que a literatura já cogitava mas nunca tinha isolado diretamente.

## Refinamento científico #3 (2026-09-05/06) — validação da extração bruta + extensão da série

Retomando os 3 pontos de fortalecimento do paper identificados na revisão crítica anterior, atacados nesta ordem: (3) validar a extração bruta → (2) estender a série além do SOL 2482 → (1) regressão direta contra pressão do REMS.

**(3) Validação da extração bruta em Python — CONCLUÍDA, pipeline original confirmado correto.**
Os arquivos RDR brutos (`.TXT`) originais usados por `mslrad_v1.py` não existem mais localmente (só sobraram os `results_detector_B/E.txt` já processados). Em vez de confiar cegamente neles, foi feita uma validação independente:
- Localizado o archive público do RAD no PDS-PPI: `https://pds-ppi.igpp.ucla.edu/data/MSL-M-RAD-3-RDR-V1.0/DATA/` — mesma estrutura de pastas `SOL_XXXXX_YYYYY` já usada no projeto, confirmando que são exatamente os mesmos produtos.
- Descoberto que cada arquivo `.TXT` (~65 MB/sol) tem um `.LBL` (PDS3 detached label, ~1.2 MB) com ponteiros de byte exatos (`^OBSnnn_TOT_DOSE_B_ELEMENT = ("arquivo.TXT", OFFSET <BYTES>)` + `START_BYTE`/`BYTES`) para cada valor de dose — isso permite extrair os valores via **HTTP Range request** direcionado, sem baixar o arquivo de 65 MB inteiro (confirmado que o servidor aceita multi-range em uma única requisição, resposta `multipart/byteranges`).
- Implementado em `python_reanalysis/rad_extract_lib.py`, validado em `python_reanalysis/04_validate_extraction.py`: reextração byte-exata de duas janelas (SOL 00000_00089 completa — 78 sols — e amostra de SOL 02359_02482 — 6 sols) comparada linha a linha contra os `results_detector_B/E.txt` locais.
- **Resultado: 3763/3763 valores comparados batem exatamente (0 discrepâncias)**, cobrindo início e fim da série original. A heurística frágil do `mslrad_v1.py` original ("valor está 2 linhas depois do header de texto") produziu dados corretos — confirmado por um método totalmente independente e byte-preciso. Gap fechado.

**(2) Extensão da série além do SOL 2482 — em andamento.**
O mesmo archive PDS-PPI tem dados do RAD disponíveis até **SOL 04710_04843** (lançado em 2026-07-29) — mais que o dobro do corte usado no paper (SOL 2482, ~set/2019). Isso cobre o mínimo solar completo (dez/2019-out/2020) e todo o ciclo solar 25 (máximo ~2024-2025) além do início do declínio, superando em cobertura temporal o próprio Ehresmann et al. (2023). Baixar os `.TXT` brutos de todos os ~2360 sols novos seria inviável (~150 GB); a extração por byte-offset via `.LBL` (mesma técnica validada acima) resolve isso, baixando só ~1.2 MB de label + poucos KB de dados por sol (~2.8 GB total). Pipeline de extração em massa sendo executado.

**(1) Regressão direta contra pressão do REMS — planejada, ainda não iniciada.**
Fonte localizada: `https://atmos.nmsu.edu/PDS/data/mslrem_1001/DATA/`, mesma estrutura de pastas `SOL_XXXXX_YYYYY`, cobertura também até SOL_04710_04843. Cada sol tem 4 produtos (`ADR`/`RMD`/`RNV`/`RTL`, 5-28 MB cada) — ainda não determinado qual contém a coluna de pressão em formato PDS4 (`.xml` + `.TAB` de registro fixo, não PDS3 com ponteiros por elemento como o RAD). Próximo passo: inspecionar o label de um desses produtos para confirmar a coluna de pressão e desenhar a mesma estratégia de Range request por linha de registro fixo.

## Refinamento científico #3 (continuação, 2026-09-06) — resultados da reanálise estendida

Os 3 pontos completados. Pipeline final: `07_backfill_original_history.py` (corrigido para lidar com sols com múltiplos arquivos de downlink — ver nota de bug abaixo) + `05_extend_series.py` → `08_consolidate.py` → `09_extended_analysis.py` + `10_solar_pressure_regression.py`. Dados brutos: `mslrad_master_sol_series.csv` (sol 0-4843, dose B/E por sol, pressão REMS por sol, número de manchas solares SILSO observado).

**Bug encontrado e corrigido durante o backfill**: 4 das 22 janelas originais (SOL_00180_00269, SOL_00939_01062, SOL_01294_01417, SOL_01649_01772) têm **múltiplos arquivos RDR para o mesmo sol** (passes de downlink diferentes — até 17 arquivos para um único sol em SOL_01649_01772). A extração inicial usava `(sol, obs_index, detector)` como chave, o que descartava observações reais de arquivos extras do mesmo sol como se fossem duplicatas. Corrigido em `rad_extract_lib.process_folder_to_csv`: `obs_index` agora é único por slot-de-arquivo dentro do sol, e o resumo de progresso é rastreado por arquivo-fonte (manifest), não por sol. Reprocessado do zero para as 4 janelas afetadas; verificado 0 chaves duplicadas após a correção.

**Bug encontrado e corrigido na consolidação**: a interpolação do número de manchas solares (SILSO) para a data de cada sol devolvia o mesmo valor (76,0 — o último ponto da série) para TODOS os 4844 sols. Causa: `pd.to_datetime(dict(...))` no pandas 2.x devolve `datetime64[us]`, enquanto `Timestamp + Timedelta` fica em `datetime64[ns]` — misturar os dois no `.astype("int64")` desalinha as unidades em 1000x, fazendo o `np.interp` saturar no último valor. Corrigido forçando `datetime64[ns]` dos dois lados antes de converter. Verificado: SSN agora mostra a variação real esperada (mínimo solar profundo perto do sol 2482 / 2019, consistente com o mínimo real do ciclo 24/25; subindo de novo até sol 4843 / mar/2026).

**Resultados (sol 1-4843, 7,24 anos marcianos — contra 3,53 no paper original):**

1. **Decomposição clássica + wavelet, série completa**: razão sazonal/residual = 0,45 (B) e 0,44 (E) — o componente sazonal continua bem menor que o ruído residual, agora com o dobro de ciclos marcianos observados. Fração da banda anual do wavelet que excede significância de 95%: 3,44% (B), 0,00% (E) — essencialmente nenhuma periodicidade anual robusta sobrevive.

2. **Regressão contínua dose ~ SSN suavizado** (substitui o teste-t binário Period1×Period2): R² = 0,64 (B) e 0,70 (E), p ≈ 0 — a atividade solar por si só explica 64-70% da variância da dose. Muito mais rigoroso que comparar duas médias de período.

3. **Controle positivo com pressão real do REMS**: a MESMA decomposição aplicada à pressão (que tem um ciclo sazonal real e bem documentado, CO2 condensando/sublimando nas calotas) dá razão sazonal/residual = 8,36 — **19x maior** que a do resíduo da dose (0,43-0,44). Isso demonstra diretamente que o método detecta sazonalidade real quando ela existe, e que a ausência de sazonalidade na dose não é um ponto cego do método.

4. **Correlação cruzada resíduo-da-dose × pressão real**: |r| máximo de apenas 0,056-0,090 (em qualquer lag entre -400 e +400 sols), sem lag fisicamente interpretável — reforça a ausência de ligação sazonal.

**Conclusão**: o resultado negativo do paper original (sazonalidade não resolvida) se mantém e fica muito mais forte — agora com o dobro dos ciclos, um teste de regressão solar contínuo e rigoroso, e um controle positivo direto com dado real de pressão mostrando que o método funciona. Isso muda a avaliação de adequação editorial (ver resposta ao usuário na conversa).

**Terceiro bug encontrado e corrigido (figuras)**: a decomposição clássica (painel "Seasonal" da Figure 6) não é robusta a SPEs — um único sol com excursão extrema cai numa fase específica da média módulo-669, e como a reconstrução tile-a esse valor médio por fase em TODOS os ciclos (`seasonal = seasonal_avg[phase]`), um evento real e pontual (ex. a SPE de maio/2024, SOL 4190 — a maior alta de radiação já registrada pelo RAD desde o pouso, segundo NASA/JPL) aparecia como um pico "sazonal" falso repetido a cada ~669 sóis por toda a figura. Corrigido com detecção robusta de outlier (mediana móvel de 21 sóis + limiar de 8 MAD) que marca e interpola sobre sóis afetados por SPE **só para o cálculo da decomposição** — o painel "Observed" e o wavelet (que já localiza esses eventos corretamente como transientes de período curto, bem abaixo da banda anual) continuam mostrando os dados brutos sem alteração. 35 sóis (B) / 29 sóis (E) marcados, agrupados em ~9 eventos distintos (consistente com os SPEs já conhecidos/citados no Table 1 original: "2014 SPE" ~sol 737-739, "2017 SPE" ~sol 1812-1816, mais eventos novos no ciclo solar 25 em 2022-2024). Figuras finais em `paper/figuras/`: `Figure6_ClassicalDecomposition_TrendSeasonalResidual_E_extended.png`, `Figure7_Wavelet_TOTALSERIES_E_extended.png`, `Figure8_SolarPressureRegression_E.png` (nova). Antigas Figure6/7 (SOL 0-2358) removidas.

**Sensibilidade da regressão solar ao mesmo evento**: a mesma SPE de maio/2024 também distorcia a regressão OLS dose~SSN suavizado (script 10) — incluí-la reduz R² de 0,84/0,85 para 0,64/0,70 (B/E) sem mudar a inclinação (~1%). Fit primário reportado exclui esse sol; ambos os números aparecem no texto/Table 4 por transparência.

**Números finais (pós-correção robusta a SPE, os que entram no paper):**
- Registro estendido: SOL 1-4843, ~13,5 anos terrestres, 7,24 anos marcianos (era 3,53).
- Validação byte-exata: 3763/3763 valores batem.
- Decomposição clássica (robusta a SPE): razão sazonal/resíduo B=0,44, E=0,37.
- Wavelet: fração da banda anual com significância >95%: B=3,44%, E=0,00%.
- OLS dose~SSN suavizado, fit primário (exclui SOL 4190): B slope=-0,0319 R²=0,840; E slope=-0,0430 R²=0,855 (p≈0 nos dois).
- OLS incluindo SOL 4190: B R²=0,639; E R²=0,702.
- Controle positivo (pressão, robusto a outlier): razão sazonal/resíduo = 9,69.
- Resíduo da dose (pós-fit solar, robusto a SPE): razão sazonal/resíduo B=0,375, E=0,317.
- Pressão mostra 25,9x (B) / 30,6x (E) mais estrutura sazonal que o resíduo da dose.
- Correlação cruzada resíduo×pressão: |r| máx B=0,089 (lag -390 sóis), E=0,146 (lag -366 sóis) — sem lag fisicamente interpretável.
- SPE de maio/2024 (SOL 4190 ≈ 2024-05-14): documentada pela NASA/JPL como a maior alta de radiação já registrada pelo RAD desde o pouso; média por sol B=71,7, E=78,4 µGy/h contra média do registro inteiro ~9,8-9,9 µGy/h.

## Atualização do texto e figuras do paper (2026-09-06)

Todo o `paper_msl_2026.docx` atualizado para refletir a série estendida (via python-docx, scripts temporários já removidos após aplicação e verificação):
- **Abstract e Introduction** reescritos com os novos números (SOL 4843, validação byte-exata, R² da regressão solar, controle positivo de pressão).
- **Métodos (2.3)**: 4 parágrafos novos (validação, extensão, regressão solar+pressão, robustez a SPE) inseridos após o parágrafo do teste-t original. Parágrafos originais sobre decomposição/wavelet (que descreviam a metodologia baseada em índice de observação do registro de 7 anos) reconciliados com uma frase apontando que as Figuras 6-7 agora usam a versão estendida com eixo de sol real.
- **Resultados 3.3**: texto/legendas atualizados para SOL 00001-04843; Table 3 recalculada na série completa; Figures 6 e 7 substituídas pelas versões estendidas (antigas removidas de `figuras/`).
- **Nova Seção 3.4** ("Solar Activity and Atmospheric Pressure Regression"): 4 parágrafos + Table 4 (nova) + Figure 8 (nova), inseridos entre 3.3 e a Discussion.
- **Discussion e Conclusion** reescritos — a lista de "future work" da Conclusion original (estender o registro, validar detector B, adicionar figura de pressão do REMS) foi cumprida nesta sessão; substituída por itens novos (modelar Φ solar diretamente, caracterizar o evento SOL 4190).
- **Referências**: adicionada SILSO World Data Center (2026) — fonte do número de manchas solares.

Verificação final: nenhuma menção órfã a "SOL 02482 (~3,5 anos marcianos)" ou "TOTALSERIES" (7 anos) sobrou fora dos dois parágrafos de reconciliação (29-30), que citam o registro original de propósito, como contexto histórico do método.

## Limpeza de auto-referência + remoção de Table 2/Figure 5 (2026-09-06)

Usuário apontou, corretamente, que o texto não deveria mencionar "análises anteriores que fizemos e não usamos mais" — isso não faz sentido num paper (que descreve o método final, não o histórico de revisão). Varredura completa:
- **Removido todo o vocabulário de auto-referência** ("earlier pass of this analysis", "original analysis", "correcting an earlier pass", "supersedes", "now extending/covering/over... (comparado a antes)") do Abstract, Introduction, Métodos 2.3, Resultados 3.3/3.4, Discussion e Conclusion. Mantidas as comparações legítimas com literatura EXTERNA (Guo, Ehresmann, Rafkin) — isso são citações reais, não auto-referência.
- **Métodos 2.3 reescrito do zero**: 7 parágrafos limpos (Table 1, extração+validação+composição do registro, decomposição clássica, wavelet, robustez a SPE, regressão solar, regressão de pressão), descrevendo só o método final, sem narrar a migração R→Python nem a extensão SOL 2482→4843 como "correção" de algo anterior.
- **Table 2 (teste-t Period1×Period2) e Figure 5 (gráfico dessa mesma comparação) removidas** — decisão confirmada com o usuário via pergunta direta: ficaram estritamente superadas pela regressão contínua dose~atividade solar (Table 3/Figure 7 novas, mesma pergunta, resposta mais rigorosa). Table 1 (médias descritivas por janela) mantida — é dado bruto de referência, não uma "análise superada". Parágrafo do teste-t em Métodos 2.3 removido junto.
- **Renumeração de Tables/Figures** após a remoção: Table 3→2, Table 4→3; Figure 6→5, 7→6, 8→7 — script de regex cobriu quase tudo automaticamente, exceto um caso com "Figures 6 and 8" (formato "and", não intervalo com traço) que passou batido e foi corrigido manualmente para "Figures 5 and 7".
- Verificação final: 0 ocorrências de qualquer termo de auto-referência remanescente; numeração de Tables (1-3) e Figures (1-7) sequencial e sem lacunas; 113 parágrafos, 3 tabelas, 7 imagens.

## Upgrade "alto impacto" (2026-09-06) — histerese, proxy físico, Lomb-Scargle, dose de missão

Usuário pediu avaliação crítica de onde o paper podia virar contribuição de alto impacto, não só publicável. Executados, na ordem sugerida:

1. **Correção de significância (autocorrelação)**: a regressão dose~SSN reincidia no mesmo problema que o próprio paper já reconhecia (teste-t não pode rodar sobre observação bruta autocorrelacionada) — agora com erro-padrão Newey-West (HAC, maxlags=400). SE inflaciona ~11x, mas mesmo assim continua extremamente significativo.
2. **Proxy físico real de GCR**: contagem do monitor de nêutrons de Oulu (NMDB, 2012-2026, baixado via query URL do nest tool) prediz a dose muito melhor que o SSN (R²~0,95-0,96 contra ~0,88).
3. **Teste de lag/histerese**: sem lag detectável entre Marte e Terra (resposta essencialmente instantânea, mesmo com suavização fina de 27 sóis). Mas **histerese real e estatisticamente significativa**: inclinação ~20% mais forte na subida do ciclo 25 que na descida do ciclo 24 (teste de interação formal, p<0,002 nos dois detectores) — comparação inédita na literatura publicada do RAD, só possível porque agora temos uma reversão completa de ciclo solar no registro.
4. **Lomb-Scargle** (terceiro método independente): potência na frequência anual não distinguível do ruído (FAP~0,99-1,00 nos dois detectores) — confirma o resultado negativo por um terceiro caminho metodológico.
5. **Dose de missão**: usando o próprio fator de qualidade do RAD (Q=3,05, Hassler et al. 2014) e a própria dose de cruzeiro deles (662 mSv), a dose de missão calculada pra cenário solar-máximo bate quase exato com o valor publicado por Hassler (983 vs 982 mSv) — validação cruzada forte do pipeline inteiro. Todos os 3 cenários (max solar, mínimo solar, ciclo 25 atual) excedem o limite de carreira da NASA (600 mSv), por 49-89%.
6. **Data/Code Availability** adicionada antes das Referências.

**5 referências novas verificadas em fonte primária** (nenhuma inventada): Scargle (1982, ApJ 263:835), Newey & West (1987, Econometrica 55:703), Mavromichalaki et al. (2011, Adv. Space Res. 47:2210, NMDB), National Academies (2021, doi:10.17226/26155, limite de 600 mSv), VanderPlas (2018, ApJS 236:16, Lomb-Scargle).

**3 novas seções** (3.5, 3.6, 3.7), **2 tabelas novas** (Table 4: regressão NM+histerese; Table 5: dose de missão), **2 figuras novas** (Figure 8: lag scan+histerese; Figure 9: periodograma Lomb-Scargle). Abstract, Discussion, Conclusion reescritos pra refletir os achados novos. Scripts em `MSL/2026/python_reanalysis/` (11_nm_lag_hysteresis.py, 12_lombscargle.py, dados em `nm_data/oulu_daily.csv`).

Estado final: 155 parágrafos, 5 tabelas, 9 imagens.

## Repositório público criado (2026-09-06)

Usuário pediu repositório git novo ("msl_2026") só com os arquivos desta última versão do manuscrito — nada das pastas antigas (`MSL/2023/`, `MSL/2026/MSLRAD BIBLIOS/` etc.). `gh` CLI não estava instalado no ambiente; baixado o binário portátil (zip da release oficial `cli/cli`, o instalador MSI falhou por exigir admin) e extraído para `%LOCALAPPDATA%\GitHubCLI`, adicionado ao PATH do usuário. Autenticado via device-code flow (usuário completou manualmente).

Duas decisões confirmadas com o usuário antes de commitar:
- **`paper/biblios/` (38MB de PDFs com copyright) excluído** do repo — risco de redistribuição mesmo em repo privado.
- **Sem rodapé "Co-Authored-By: Claude"** nos commits — usuário confirmou que a preferência histórica (nunca incluir) vale também aqui, apesar da configuração desta sessão pedir o contrário.

Repositório: **https://github.com/cesarrlamaral/msl_2026** (privado). Estrutura: `manuscript/` (docx + EDITORIAL_NOTES.md + figuras/) e `code/` (scripts 01-12 + rad_extract_lib.py/rems_extract_lib.py + `code/data/` com os datasets derivados pequenos — série mestre por sol, contagem do Oulu, CSVs de resumo). Excluídos do repo: PDFs de bibliografia, exports brutos por observação (extended_series/, pressure_series/, decomposition_annual_period.csv — muito grandes, a série mestre já agrega o essencial), logs e `__pycache__`. O Data and Code Availability Statement do paper foi atualizado do texto genérico "disponível mediante solicitação" pra apontar direto pro link real do repositório.

## Formatação completa para submissão ao JGR: Planets / AGU (2026-09-06)

Usuário pediu formatação completa (texto, figuras, tabelas, referências, autoria) seguindo as normas da AGU (editora do JGR: Planets). Pesquisa feita direto no site da AGU (agu.org bloqueava fetch automatizado por proteção anti-bot — contornado usando o navegador Claude-in-Chrome) nas páginas oficiais "Grammar and Style Guide" e "Text & Graphics Requirements".

**Mudanças aplicadas:**
- **Autoria/afiliação**: autores com sobrescrito numérico (Cesar Amaral¹,*, Dafne Adriana Abreu dos Anjos¹, Anna Luisa dos Santos Donato¹, Leticia Bastos Eller¹), afiliação completa da UERJ, linha de autor correspondente com e-mail.
- **Key Points**: 3 pontos, todos ≤140 caracteres (obrigatório AGU).
- **Abstract**: cortado de 420 para 249 palavras (limite AGU é <250), citações removidas (regra AGU: evitar citação no abstract a menos que essencial).
- **Plain Language Summary**: novo, obrigatório pro JGR:Planets especificamente, 200 palavras, sem jargão/siglas.
- **Keywords**: 6 palavras-chave livres adicionadas.
- **Tables e Figures movidas pro final do documento** (depois de References), ordem exigida pela AGU: Text → Acknowledgements → Open Research → References → Tables → Figures. As citações no corpo do texto ("Table 3", "Figure 5, trend panel") continuam intactas.
- **Conflict of Interest** (declaração padrão) e **Acknowledgments** (com placeholder pro usuário preencher financiamento) adicionados na ordem certa.
- **"Data and Code Availability" renomeada pra "Open Research"** (nome exigido pela AGU) e conectada a citações formais na lista de Referências (regra AGU: dado/código citado no texto precisa ter entrada na lista de Referências).
- **Referências**: todas as 28 entradas reordenadas alfabeticamente por primeiro autor (regra letra-por-letra da AGU), DOIs convertidos de `doi:10.xxxx` pra `https://doi.org/10.xxxx` (formato exigido), autores de 8+ truncados pra 6+"et al." (Ehresmann et al. 2014 tinha 14 autores — busquei a lista completa no PDF original em `paper/biblios/` pra truncar corretamente; Guo et al. 2021 tinha 10), títulos convertidos pra sentence case, nome de periódico+volume em itálico. Adicionadas 4 referências novas de dado/software (NASA PDS RAD, NASA PDS REMS, SILSO dataset, código no GitHub) exigidas pela seção Open Research.
- **Espaçamento duplo** (regra AGU) aplicado no estilo Normal + numeração contínua de linha (`w:lnNumType`) na seção do documento — ambos exigidos pela submissão.

**Bug corrigido no meio do caminho**: minhas primeiras tentativas de mover Tables/Figures pro final usaram `body.append()` bruto, que colocou o conteúdo DEPOIS do elemento `w:sectPr` (que precisa ser sempre o último filho do body) — documento tecnicamente inválido. Corrigido reconstruindo com `sectPr.addprevious()`.

**Pendências que só o usuário pode resolver** (documentadas na seção "Ainda em aberto" abaixo).

## Ainda em aberto

- **Comparação quantitativa com o campo magnético terrestre.** O rascunho original tinha um "XX%" não preenchido comparando o campo de Marte com o da Terra (Discussion). Removido do corpo do texto por não fazer sentido como placeholder numa versão "pronta", mas o ponto de fundo é real: Marte não tem campo global comparável (só magnetismo remanescente crustal localizado, ex. Terra Cimmeria ~1600 nT, contra ~25-65 µT do campo superficial terrestre) — uma comparação percentual direta não é fisicamente bem definida do jeito que estava formulada. Se quiser esse dado no texto, precisa ser reformulado como comparação crustal-vs-global, não uma razão simples.
- **6 referências reais disponíveis em `paper/biblios/` mas não citadas no texto atual** — candidatas a incorporar se expandir Métodos/Discussion: Zeitlin et al. (2016, calibração do RAD — candidata natural pra 2.1 RAD Instrument), Matthiä & Berger (2017, GEANT4/PLANETOCOSMICS), Matthiä et al. (workshop de modelos, LSSR), Kim et al. (comparação Badhwar-O'Neill/HZETRN), Gronoff, Norman & Mertens (2015, HZETRN vs. Planetocosmics), Ratliff, Smith & Heilbronn (simulação MCNP6).

## Resolvido nesta sessão

- **ORCID, financiamento e conflito de interesse preenchidos**: ORCID do autor correspondente (https://orcid.org/0000-0002-4314-3517) adicionado à linha de autor correspondente na página de título; Acknowledgments atualizado com declaração explícita de que a pesquisa não recebeu financiamento externo dedicado (placeholder removido); declaração de "no conflicts of interest" confirmada como válida para os 4 autores. Único passo que ainda depende do sistema de submissão da AGU (não do manuscrito): vincular esse ORCID à conta usada no GEMS no momento da submissão.
- **Depósito arquivado no Zenodo com DOI persistente**: código + dados derivados (README, `code/`) publicados como depósito no Zenodo — https://doi.org/10.5281/zenodo.22550046 — satisfazendo a exigência da AGU de "repositório confiável". Seção Open Research e a entrada de referência do software atualizadas para citar o DOI do Zenodo (agora com os 4 autores, no lugar do link simples do GitHub, que continua mencionado apenas como material suplementar de desenvolvimento).
- **Equação da dose (Seção 2.2) reconstruída como objeto de equação nativo do Word (OMML)**, substituindo a imagem PNG (`Figure_Dose_Equation.png`, agora sem uso no documento — arquivo mantido na pasta `figuras/` apenas como registro histórico). A equação — soma dupla (j, área) + integral dupla (ε_min a ε_max, 0 a ∞) de λ_j(E,ε)F_j(Φ,P,E) dEdε/m — foi montada diretamente em XML (`m:oMath`/`m:nary`/`m:sSub`/`m:d`), com dois tab stops (centralizado + direita) reproduzindo o layout padrão do Word para equação numerada, mantendo o rótulo "Equation 1." alinhado à margem direita. Chamada no texto corrigida de "the equation 1" para "Equation 1" (sem artigo, no padrão AGU de citação de objetos numerados). Validado reabrindo o `.docx` via python-docx e checando bom formação do XML após a edição.
- **Figura de pressão do REMS descartada (não era essencial).** Avaliação: o parágrafo (Seção 2.2) só descreve contexto já publicado por Harri et al. (2014) — não é resultado original deste paper — e a Figura 1 já tem um painel (c) de pressão cobrindo a variação diurna que sustenta o argumento da Seção 3.1. Decisão: manter o parágrafo como prosa bem citada, sem exigir uma figura própria. `[FIGURE PENDING]` removido do texto.
- **Referências completas de Roesler (1998) e Zeitlin (2013)** — volume/páginas/DOI confirmados via busca e preenchidos (antes tinham nota "não confirmado nesta passada").
- **Citação "Harri et al. [2014]" (pressão do REMS)** — inicialmente sinalizada como possível erro (só um Harri 2014 sobre umidade tinha sido encontrado). Busca adicional achou que Harri publicou DOIS papers em 2014 — a citação já estava certa desde o início. Referência completa adicionada: Harri et al. (2014), "Pressure observations by the Curiosity rover: Initial results," JGR Planets 119(1), 82-92, doi:10.1002/2013JE004423.
- **Frase sobre HZETRN2015 na Discussion** — parafraseava o abstract do Guo et al. (2017) como se fosse método deste paper. Reescrita como comparação ao achado do Guo et al. (2017), não como afirmação sobre o que este paper fez.
- **3 referências não localizadas mesmo após múltiplas buscas** (Richter & Rasch 2008 — nomes possivelmente mal transcritos no rascunho original; Nakamura et al. 1987; Florek et al. 1996) — todas citadas originalmente só pra sustentar o mesmo ponto que Roesler et al. (1998)/Bazilevskaya & Svirzhevskaya (1998) já cobrem, ambas verificadas. Decisão editorial: removidas as citações não verificáveis do texto (Seção 3.1), mantidas só as duas fontes confirmadas. Nenhuma referência foi inventada para preencher a lacuna.
- **Dump cru de `summary()` do R na Seção 3.3** (usuário apontou) — virou Table 3 formatada, com precisão arredondada a 3 casas decimais.
- **Tabela 1**: vírgula decimal (padrão BR) trocada por ponto, precisão padronizada a 3 casas decimais, corrigido typo de label ("SOL 016490_01772" → "SOL 01649_01772").
- **Abstract e Conclusion escritos** a partir dos achados já estabelecidos e citados no corpo do paper (não são mais placeholder `[TODO]`).
- **Reanálise científica completa em Python** (decomposição clássica + wavelet com escala anual real, teste-t reprodutível) — ver detalhe abaixo.

## Refinamento científico — reanálise em Python, migração de R

R não estava instalado no ambiente; usuário pediu pra não usar mais R (ver `feedback_mslrad_python_nao_r`). Scripts em `MSL/2026/python_reanalysis/`:
- `01_seasonal_decompose.py` — decomposição clássica (tendência/sazonal/resíduo) via média móvel, com período real de 1 ano marciano (67.636 observações ≈ 668,6 sóis), corrigindo o `frequency=100` (≈1 sol) do R original. STL do statsmodels (com e sem `robust=True`) travou/não terminou em tempo viável com essa janela gigante — implementação manual via média móvel (mesmo algoritmo do `decompose()` do R) resolveu em segundos.
- `02_wavelet_totalseries.py` — wavelet de Morlet (Torrence & Compo 1998, via `pycwt`) na série completa de 7 anos, com faixa de período estendida até ~786 sóis, corrigindo o `upperPeriod=500` (≈5 sóis) do R original.
- `03_period_ttest.py` — teste-t de Welch reprodutível comparando Period 1 vs Period 2 sobre as 22 médias de janela (`medias.xlsx`), já que o teste original nunca foi salvo em nenhum script R. Resultado bate na ordem de grandeza e na conclusão (ambos p ≪ 0,05), mas não no valor exato — a variante exata usada originalmente não pôde ser confirmada.

**Achado central**: mesmo com a escala corrigida, nem o wavelet nem a decomposição clássica mostram um ciclo sazonal anual limpo na série de 7 anos (~3,5 anos marcianos) — o que aparece é o mesmo efeito de borda ("domo") característico de uma tendência secular forte, não uma banda periódica sustentada. **Conclusão adotada no paper**: o aumento de dose de longo prazo é uma tendência suave ligada à fase decrescente do ciclo solar, não um ciclo sazonal marciano resolvido — um resultado negativo defensável, incorporado ao texto (Introduction, Seção 3.3, Discussion, Conclusion).

Duas figuras novas (`Figure6_ClassicalDecomposition...png`, `Figure7_Wavelet_TOTALSERIES...png`) em `paper/figuras/`, embutidas na Seção 3.3 do docx.

## Fontes descartadas da bibliografia (não citáveis como referência formal)
Ficaram em `2026/MSLRAD BIBLIOS/` — não foram movidas para `paper/biblios/`:
- Palestra de D. Hassler (Space Weather at Mars, slides de apresentação).
- Resumo de congresso LPSC ("New Results from the MSL-RAD Experiment on Curiosity").
- Notícia da revista *Science* (jornalismo, não pesquisa original).
- Relatório técnico da NASA ("Updates from the MSL-RAD Experiment").
