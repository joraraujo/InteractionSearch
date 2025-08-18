import os
import itertools
from langchain_google_genai import ChatGoogleGenerativeAI
from Bio import Entrez
import streamlit as st
import pandas as pd

# =============================
# Configurações iniciais
# =============================
google_api_key = st.secrets["GOOGLE_API_KEY"]
Entrez.email = "jorge.oa@hotmail.com"

# Inicializa LLM
llm = ChatGoogleGenerativeAI(model="models/gemini-2.5-flash", temperature=0.2)

# Configuração da página
st.set_page_config(page_title="Pesquisa de Interações", layout="centered")
st.title("🔍 Pesquisa de Interações")

# Histórico do chat
if "messages" not in st.session_state:
    st.session_state.messages = []

# Exibe histórico
for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])

# =============================
# Funções auxiliares
# =============================

def buscar_artigos_pubmed(a1, a2, max_artigos=5):
    """
    Busca artigos no PubMed e retorna lista de dicts com id, resumo e link.
    """
    query = f"{a1} AND {a2}"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_artigos)
    record = Entrez.read(handle)
    ids = record.get("IdList", [])

    artigos = []
    for id in ids:
        try:
            handle = Entrez.efetch(db="pubmed", id=id, rettype="abstract", retmode="text")
            resumo = handle.read().strip()
        except Exception:
            resumo = "Resumo não disponível"
        artigos.append({
            "id": id,
            "resumo": resumo if resumo else "Resumo não disponível",
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{id}"
        })
    return artigos

def markdown_to_df(md_text: str) -> pd.DataFrame:
    colunas = [
        "Substâncias envolvidas",
        "Existe interação? (sim/não)",
        "Tipo de interação",
        "Forma farmacêutica",
        "Link da fonte",
    ]

    linhas_validas = []
    for raw in md_text.splitlines():
        line = raw.strip()
        if not line or "|" not in line:
            continue
        if "Substâncias envolvidas" in line or set(line.replace("|", "").strip()) == {"-"}:
            continue

        cells = [c.strip() for c in line.strip("|").split("|")]
        cells = [c for c in cells if c != ""]  # remove células vazias

        if len(cells) < 5:
            cells += [""] * (5 - len(cells))
        elif len(cells) > 5:
            cells = cells[:4] + [" | ".join(cells[4:])]

        # correção de deslocamento
        if cells[0] == "" and len(cells) > 1 and " + " in cells[1]:
            cells = [cells[1]] + cells[2:] + [""]

        linhas_validas.append(cells[:5])

    df = pd.DataFrame(linhas_validas, columns=colunas)
    if df.empty:
        return df
    df = df.replace({"": "não informado"})
    return df

def prompt_por_artigo(par, artigo):
    """
    Monta prompt para análise de um único artigo.
    """
    a1, a2 = par
    prompt = (
        "Você é um especialista em química/farmacêutica. ANALISE APENAS o resumo abaixo.\n\n"
        f"Par de substâncias: {a1} + {a2}\n"
        f"Link do artigo: {artigo['link']}\n\n"
        "Resumo do artigo:\n"
        f"{artigo['resumo']}\n\n"
        "INSTRUÇÕES DE SAÍDA (muito importantes):\n"
        "- Responda com exatamente UMA LINHA no formato Markdown, usando '|' como separador,\n"
        "  sem '|' no início nem no fim da linha.\n"
        "- A linha deve conter 5 campos (nesta ordem):\n"
        "  Substâncias envolvidas | Existe interação? (sim/não) | Tipo de interação | Forma farmacêutica | Link da fonte\n"
        "- Para 'Substâncias envolvidas' use: '<A1> + <A2>' (por exemplo: Ácido ascórbico + Riboflavina).\n"
        "- Para 'Link da fonte' use o link exato fornecido acima.\n"
        "- Se não houver interação descrita, não liste na tabela.\n"
        "- Se algum campo não puder ser inferido do artigo, escreva 'não informado' nesse campo.\n"
        "- NÃO inclua cabeçalho, linhas de separador (---) nem texto adicional.\n"
    )
    return prompt

def gerar_tabela_interacoes(ativos, max_por_par=5):
    """
    Para cada par: busca os artigos e pede ao LLM que analise cada artigo separadamente.
    Retorna DataFrame consolidado.
    """
    pares = list(itertools.combinations(ativos, 2))
    df_total = pd.DataFrame(columns=[
        "Substâncias envolvidas",
        "Existe interação? (sim/não)",
        "Tipo de interação",
        "Forma farmacêutica",
        "Link da fonte",
    ])

    total_pairs = len(pares)
    pair_idx = 0

    for par in pares:
        pair_idx += 1
        a1, a2 = par
        st.write(f"🔎 Buscando artigos para par ({pair_idx}/{total_pairs}): **{a1} + {a2}**")
        artigos = buscar_artigos_pubmed(a1, a2, max_artigos=max_por_par)

        if not artigos:
            df_total = pd.concat([df_total, pd.DataFrame([{
                "Substâncias envolvidas": f"{a1} + {a2}",
                "Existe interação? (sim/não)": "não informado",
                "Tipo de interação": "não informado",
                "Forma farmacêutica": "não informado",
                "Link da fonte": "N/A",
            }])], ignore_index=True)
            continue

        artigo_idx = 0
        for artigo in artigos:
            artigo_idx += 1
            st.write(f" • [{artigo_idx}/{len(artigos)}] analisando: {artigo['link']}")
            prompt = prompt_por_artigo(par, artigo)
            try:
                resposta = llm.invoke(prompt).content
            except Exception as e:
                resposta = ""

            df_parsed = markdown_to_df(resposta)
            if not df_parsed.empty:
                df_parsed.loc[:, "Link da fonte"] = artigo['link']  # Garante link correto
                df_total = pd.concat([df_total, df_parsed], ignore_index=True)
            else:
                df_total = pd.concat([df_total, pd.DataFrame([{
                    "Substâncias envolvidas": f"{a1} + {a2}",
                    "Existe interação? (sim/não)": "não informado",
                    "Tipo de interação": "não informado",
                    "Forma farmacêutica": "não informado",
                    "Link da fonte": artigo['link'],
                }])], ignore_index=True)

    df_total = df_total.fillna("não informado")
    return df_total

# =============================
# UI Streamlit
# =============================

with st.sidebar:
    st.header("Configurações")
    max_por_par = st.slider("Máx. artigos por par", min_value=1, max_value=20, value=5, step=1,
                            help="Número máximo de artigos PubMed a buscar por cada par de ativos.")
    temperature = st.slider("Temperatura do LLM", min_value=0.0, max_value=1.0, value=0.2, step=0.1)
    llm.temperature = temperature

st.markdown("Digite os ativos separados por vírgula no campo de chat abaixo.")

# Entrada de ativos (chat)
if prompt := st.chat_input("Informe os ativos separados por vírgula..."):
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    ativos = [a.strip() for a in prompt.split(",") if a.strip()]

    if len(ativos) < 2:
        resposta = "⚠️ Informe pelo menos dois ativos."
        with st.chat_message("assistant"):
            st.markdown(resposta)
        st.session_state.messages.append({"role": "assistant", "content": resposta})
    else:
        with st.spinner("Executando buscas e análises (isso pode levar um tempo)..."):
            df = gerar_tabela_interacoes(ativos, max_por_par=max_por_par)

        with st.chat_message("assistant"):
            st.dataframe(df, use_container_width=True)

        st.session_state.messages.append({"role": "assistant", "content": "Segue a tabela (uma linha por artigo)."})
