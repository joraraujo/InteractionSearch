import os
import itertools
from langchain_google_genai import ChatGoogleGenerativeAI
from Bio import Entrez
import streamlit as st
import pandas as pd

# =============================
# Configurações iniciais
# =============================
os.environ["GOOGLE_API_KEY"] = "AIzaSyCNQbHda-5qjwsKnkVK_S7N9aW8jRZ488c"
Entrez.email = "jorge.oa@hotmail.com"

# Inicializa LLM
llm = ChatGoogleGenerativeAI(model="models/gemini-2.5-flash", temperature=0.3)

# Configuração da página
st.set_page_config(page_title="Pesquisa de Interações", layout="centered")
st.title("🔍 Pesquisa de Interações entre Ativos")

# =============================
# Histórico do chat
# =============================
if "messages" not in st.session_state:
    st.session_state.messages = []

# Exibe histórico
for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])

# =============================
# Funções auxiliares
# =============================

def buscar_artigos_pubmed(a1, a2, max_artigos=3):
    """Busca artigos no PubMed e retorna abstracts + links."""
    query = f"{a1} AND {a2}"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_artigos)
    record = Entrez.read(handle)
    ids = record["IdList"]

    artigos = []
    for id in ids:
        handle = Entrez.efetch(db="pubmed", id=id, rettype="abstract", retmode="text")
        resumo = handle.read().strip()
        artigos.append({
            "id": id,
            "resumo": resumo if resumo else "Resumo não disponível",
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{id}"
        })
    return artigos

def gerar_tabela_interacoes(ativos):
    """Busca abstracts no PubMed e gera análise com o Gemini."""
    pares = list(itertools.combinations(ativos, 2))
    prompt_llm = (
        "Você é um especialista em química farmacêutica.\n"
        "Com base nos resumos dos artigos a seguir, identifique se há interação entre os pares de substâncias listados.\n"
        "Responda estritamente no formato de tabela Markdown com as colunas:\n"
        "Substâncias envolvidas | Existe interação? (sim/não) | Tipo de interação "
        "(química, física, fotossensibilidade, etc.) | Forma farmacêutica | Link da fonte\n\n"
    )

    for a1, a2 in pares:
        artigos = buscar_artigos_pubmed(a1, a2)
        if artigos:
            prompt_llm += f"\n### Par: {a1} + {a2}\n"
            for artigo in artigos:
                prompt_llm += f"- Link: {artigo['link']}\nResumo: {artigo['resumo']}\n"
        else:
            prompt_llm += f"\n### Par: {a1} + {a2}\nNenhum artigo encontrado.\n"

    # Chamada ao LLM
    resposta_llm = llm.invoke(prompt_llm).content

    # Converter resposta em DataFrame
    linhas = [l.strip() for l in resposta_llm.splitlines() if l.strip() and "|" in l]
    dados = []
    for linha in linhas:
        partes = [p.strip() for p in linha.split("|")]
        while len(partes) < 5:  # garantir sempre 5 colunas
            partes.append("")
        dados.append(partes[:5])

    colunas = ["Substâncias envolvidas", "Existe interação? (sim/não)",
               "Tipo de interação", "Forma farmacêutica", "Link da fonte"]
    return pd.DataFrame(dados, columns=colunas)

# =============================
# Entrada de ativos (chat)
# =============================
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
        df = gerar_tabela_interacoes(ativos)

        # Exibe no chat
        with st.chat_message("assistant"):
            st.dataframe(df, use_container_width=True)

        st.session_state.messages.append({"role": "assistant", "content": "Segue a tabela com as interações encontradas."})
