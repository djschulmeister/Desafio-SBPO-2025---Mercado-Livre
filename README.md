# Desafio-SBPO-2025---Mercado-Livre
Alguns detalhes para a execução das instâncias

Nossa equipe não conseguiu executar o arquivo "run_challenge", então as instâncias foram executadas manualmente via GitBash, utilizando os seguintes comandos:

cd /c/[CAMINHO_DO_PROJETO]/challenge-sbpo-2025

mvn clean package  

# Execução das instâncias:  
java -jar target/ChallengeSBPO2025-1.0.jar datasets/a/Instance_0001.txt Resultados/solution1.txt  
java -jar target/ChallengeSBPO2025-1.0.jar datasets/a/Instance_0002.txt Resultados/solution2.txt  
...  

# Validação com checker.py:  
python checker.py datasets/a/Instance_0001.txt Resultados/solution1.txt  
python checker.py datasets/a/Instance_0002.txt Resultados/solution2.txt  
... 

Para organização, criamos uma pasta chamada Resultados para armazenar todos os arquivos de solução (solution1.txt, solution2.txt, etc.). 
Além disso, geramos um arquivo RUN&TEST.txt dentro dessa pasta com os comandos utilizados, para garantir rastreabilidade.
