# CPLEX Python isolado

Ambiente: `.venv-cplex/` na raiz do projeto, Python3.10arm64.
Instalação oficial via PyPI: `cplex` e `docplex`.

```sh
source .venv-cplex/bin/activate
python experiments/cplex_installation/check_installation.py
```

A verificação resolve o mesmo pequeno MIP pelas APIs CPLEX e DOcplex, confere
ótimo3 e testa se o runtime rejeita um modelo trivial com1001variáveis.
O resultado real fica em installation_check.json. Versões exatas devem ser
lidas desse arquivo e do requirements instalado, não presumidas.

Este ambiente não altera as dependências nem o solver dos experimentos PCCTP
em execução. Instalar o runtime não altera automaticamente o código C++/Gurobi
do repositório. O pacote Python não inclui a IDE nem o SDK C++ do Optimization
Studio. Não chamar este ambiente de edição acadêmica/completa sem ativação.

Documentação oficial de instalação e atualização para uma instalação licenciada:
https://www.ibm.com/docs/en/icos/22.1.2?topic=rnco2-cplex-python-api-now-installed-using-pip-conda
A edição Community impõe limites de tamanho. Para modelos maiores, utilizar
uma instalação/licença IBM acadêmica ou comercial apropriada; não contornar
os limites com remendos ou declarar que todos os benchmarks podem ser resolvidos.
