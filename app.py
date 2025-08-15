import os
import itertools
from langchain_google_genai import ChatGoogleGenerativeAI
from Bio import Entrez
import streamlit as st
import pandas as pd
import io

# Configurações iniciais
os.environ["GOOGLE_API_KEY"] = "AIzaSyCNQbHda-5qjwsKnkVK_S7N9aW8jRZ488c"
Entrez.email = "jorge.oa@hotmail.com"

# Inicializa LLM
llm = ChatGoogleGenerativeAI(model="models/gemini-2.5-flash", temperature=0.3)

# Configuração da página
st.set_page_config(page_title="Pesquisa de Interações", layout="centered")
st.title("🔍 Pesquisa de Interações por Pares de Ativos")

# Histórico do chat
if "messages" not in st.session_state:
    st.session_state.messages = []

# Exibe histórico
for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])

# Entrada de ativos
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
        pares = list(itertools.combinations(ativos, 2))
        pares_com_publicacoes = []

        for a1, a2 in pares:
            query = f"{a1} {a2} interaction"
            handle = Entrez.esearch(db="pubmed", term=query)
            record = Entrez.read(handle)
            count = int(record["Count"])

            if count > 0:
                link = f"https://pubmed.ncbi.nlm.nih.gov/?term={a1.replace(' ', '+')}+{a2.replace(' ', '+')}"
                pares_com_publicacoes.append((a1, a2, link))
            else:
                pares_com_publicacoes.append((a1, a2, ""))

        # Monta prompt único para o Gemini
        prompt_llm = (
            "Analise os seguintes pares de substâncias e informe se existe interação entre elas. "
            "Responda estritamente no formato de tabela Markdown com as colunas: "
            "Substâncias envolvidas | Existe interação? (sim/não) | Tipo de interação "
            "(química, física, fotossensibilidade, etc.) | Forma farmacêutica | Link da fonte\n\n"
        )

        for a1, a2, link in pares_com_publicacoes:
            prompt_llm += f"- {a1} + {a2} | Link: {link if link else 'N/A'}\n"

        # Consulta única ao LLM
        resposta_llm = llm.invoke(prompt_llm).content

        # Limpeza automática da resposta e conversão para DataFrame
        linhas = [l.strip() for l in resposta_llm.splitlines() if l.strip() and "|" in l]
        dados = []
        for linha in linhas:
            partes = [p.strip() for p in linha.split("|")]
            # Garantir que sempre tenha 5 colunas
            while len(partes) < 5:
                partes.append("")
            dados.append(partes[:5])

        colunas = ["Substâncias envolvidas", "Existe interação? (sim/não)", "Tipo de interação",
                   "Forma farmacêutica", "Link da fonte"]
        df = pd.DataFrame(dados, columns=colunas)

        # Exibe no chat
        with st.chat_message("assistant"):
            st.dataframe(df, use_container_width=True)

        st.session_state.messages.append({"role": "assistant", "content": "Segue a tabela com as interações encontradas."})
