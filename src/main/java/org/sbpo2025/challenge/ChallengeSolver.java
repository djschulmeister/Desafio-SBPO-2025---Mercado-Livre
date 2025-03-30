package org.sbpo2025.challenge;

import org.apache.commons.lang3.time.StopWatch;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.HashMap;
import java.util.Collections;
import java.util.Comparator;
import ilog.concert.*;
import ilog.cplex.*;






public class ChallengeSolver {
    private final long MAX_RUNTIME = 600000; // milliseconds; 10 minutes

    protected List<Map<Integer, Integer>> orders;
    protected List<Map<Integer, Integer>> aisles;
    protected int nItems;
    protected int waveSizeLB;
    protected int waveSizeUB;

    public ChallengeSolver(
            List<Map<Integer, Integer>> orders, List<Map<Integer, Integer>> aisles, int nItems, int waveSizeLB, int waveSizeUB) {
        this.orders = orders;
        this.aisles = aisles;
        this.nItems = nItems;
        this.waveSizeLB = waveSizeLB;
        this.waveSizeUB = waveSizeUB;
    }


public static class OptimizationResult {
    private final double objectiveValue;
    private final List<Integer> selectedOrders;
    private final List<Integer> selectedAisles;

    public OptimizationResult(double objectiveValue, 
                            List<Integer> selectedOrders, 
                            List<Integer> selectedAisles) {
        this.objectiveValue = objectiveValue;
        this.selectedOrders = selectedOrders;
        this.selectedAisles = selectedAisles;
    }

    public double getObjectiveValue() { return objectiveValue; }
    public List<Integer> getSelectedOrders() { return selectedOrders; }
    public List<Integer> getSelectedAisles() { return selectedAisles; }
}

public static OptimizationResult runModel(
    double CB, int LB, int UB,
    int NO, int NA, int NI, int L,
    int[] F,
    List<List<List<Integer>>> A,
    List<List<Integer>> Ia,
    List<List<List<Integer>>> O,
    List<List<Integer>> Io,
    List<Integer> CI,
    List<List<Integer>> G,
    List<List<Integer>> fff,
    List<Integer> ITEM_N,
    Set<Integer> O_N,
    double timeLimit) {
        
    IloCplex cplex = null;
    try {
        cplex = new IloCplex();
        
        // 1. Variáveis de decisão
        IloIntVar[] p = cplex.boolVarArray(NO);
        IloIntVar[] c = cplex.boolVarArray(NA);  // Corrigido para NA


// Definindo nomes para as variáveis de p
for (int i = 0; i < NO; i++) {
    p[i].setName("p_" + i);
}

// Definindo nomes para as variáveis de c
for (int j = 0; j < NA; j++) {
    c[j].setName("c_" + j);
}
        // 2. Função objetivo
        IloLinearNumExpr objective = cplex.linearNumExpr();
        for (int o = 0; o < NO; o++) {
            objective.addTerm((double) F[o]/L, p[o]); // Cast para double
        }
        cplex.addMaximize(objective);

        // 3. Restrições
        // 3.1 Corredores fixos
        for (int a : CI) {
            if (a < NA) cplex.addEq(c[a], 0);
        }

        // 3.2 Relação pedidos-corredores
        for (int i = 0; i < G.size(); i++) {
            List<Integer> group = G.get(i);
            if (!group.isEmpty() && i < NA) {
                IloLinearNumExpr expr = cplex.linearNumExpr();
                for (int o : group) {
                    expr.addTerm(1.0, p[o]);
                }
                cplex.addEq(expr, cplex.prod(group.size(), c[i]));
            }
        }

        // 3.3 Pedidos excluídos
        for (int o : O_N) {
            if (o < NO) cplex.addEq(p[o], 0);
        }

        // 3.4 Restrição fff
        for (int i = 0; i < fff.size(); i++) {
            List<Integer> group = fff.get(i);
            if (!group.isEmpty()) {
                IloLinearNumExpr expr = cplex.linearNumExpr();
                for (int o : group) {
                    expr.addTerm(1.0, p[o]);
                }
                cplex.addGe(expr, cplex.prod(group.size(), c[i]));
            }
        }

        // 3.5 Limites de produção
        IloLinearNumExpr productionExpr = cplex.linearNumExpr();
        for (int o = 0; o < NO; o++) {
            productionExpr.addTerm(F[o], p[o]);
        }
        cplex.addGe(productionExpr, CB * L + 1);
        cplex.addGe(productionExpr, LB);
        cplex.addLe(productionExpr, UB);  // Corrigido para addLe

       // 3.6 Capacidade dos itens (versão corrigida)
for (int i = 0; i < NI; i++) {
    if (!ITEM_N.contains(i)) {
        // Demanda total do item i nos pedidos selecionados
        IloLinearNumExpr demandaTotal = cplex.linearNumExpr();
        for (int o : Io.get(i)) {  // Pedidos que usam o item i
            int qtdPedido = getItemQuantity(O, o, i);
            demandaTotal.addTerm(qtdPedido, p[o]);
        }
        


        // Capacidade total do item i nos corredores selecionados
        IloLinearNumExpr capacidadeTotal = cplex.linearNumExpr();
        for (int a : Ia.get(i)) {  // Corredores que contêm o item i
            int capCorredor = getAisleCapacity(A, a, i);
            capacidadeTotal.addTerm(capCorredor, c[a]);
        }
        
        // Restrição: Σ(qtd_pedido * p_o) <= Σ(cap_corredor * c_a)
        cplex.addLe(demandaTotal, capacidadeTotal);
    }
}



        // 3.7 Número de corredores
        IloLinearNumExpr aisleCount = cplex.linearNumExpr();
        for (int a = 0; a < NA; a++) {
            aisleCount.addTerm(1.0, c[a]);
        }
        cplex.addEq(aisleCount, L);

        // 4. Configurações do solver
        cplex.setParam(IloCplex.Param.TimeLimit, timeLimit);
        cplex.setParam(IloCplex.Param.MIP.Display, 2);


cplex.exportModel("modelo.lp");

        // 5. Resolver
        if (cplex.solve()) {
            List<Integer> selectedOrders = new ArrayList<>();
            List<Integer> selectedAisles = new ArrayList<>();

            for (int o = 0; o < NO; o++) {
                if (cplex.getValue(p[o]) > 0.9) selectedOrders.add(o);
            }

            for (int a = 0; a < NA; a++) {
                if (cplex.getValue(c[a]) > 0.9) selectedAisles.add(a);
            }

            return new OptimizationResult(
                cplex.getObjValue(),
                selectedOrders,
                selectedAisles
            );
        }
        return new OptimizationResult(0.0, Collections.emptyList(), Collections.emptyList());
        
    } catch (IloException e) {
        System.err.println("Erro CPLEX: " + e.getMessage());
        return new OptimizationResult(0.0, Collections.emptyList(), Collections.emptyList());
    } finally {
        if (cplex != null) cplex.end();
    }
}

// Métodos auxiliares necessários
private static int getItemQuantity(List<List<List<Integer>>> orders, int orderIndex, int itemId) {
    for (List<Integer> pair : orders.get(orderIndex)) {
        if (pair.get(0) == itemId) {
            return pair.get(1);
        }
    }
    return 0;
}

private static int getAisleCapacity(List<List<List<Integer>>> aisles, int aisleIndex, int itemId) {
    for (List<Integer> pair : aisles.get(aisleIndex)) {
        if (pair.get(0) == itemId) {
            return pair.get(1);
        }
    }
    return 0;

}  // Classe interna para estado da heurística
    private static class HeuristicState {
        double FO;
        Object ySolution;
        List<Double> SolF = new ArrayList<>();
        List<Integer> SolL = new ArrayList<>();
        double LMax;
        int L;
    }



   // Classe de resultado do pré-processamento (mantida única)
    public static class PreprocessingResult {
        public List<List<List<Integer>>> A;
        public int[] d;
        public int[] CAP;
        public List<Set<Integer>> ITEM_A;
        public List<Set<Integer>> ITEM_B;
        public List<List<Integer>> DOM;
        public Set<Integer> aux;
        public List<Set<Integer>> ITEM_O;
        public List<List<Integer>> fff;
        public List<List<Integer>> ggg;
        public List<List<Integer>> G;

        public PreprocessingResult(
            List<List<List<Integer>>> A, 
            int[] d, 
            int[] CAP,
            List<Set<Integer>> ITEM_A,
            List<Set<Integer>> ITEM_B,
            List<List<Integer>> DOM,
            Set<Integer> aux,
            List<Set<Integer>> ITEM_O,
            List<List<Integer>> fff,
            List<List<Integer>> ggg,
            List<List<Integer>> G) {
            
            this.A = A;
            this.d = d;
            this.CAP = CAP;
            this.ITEM_A = ITEM_A;
            this.ITEM_B = ITEM_B;
            this.DOM = DOM;
            this.aux = aux;
            this.ITEM_O = ITEM_O;
            this.fff = fff;
            this.ggg = ggg;
            this.G = G;
        }
    }


//------------------------------------------------------------


//-------------------------------------------------------------

//------------------------------------------------------------


//-------------------------------------------------------------

//------------------------------------------------------------


//-------------------------------------------------------------


 // Método principal para resolver o problema

    public ChallengeSolution solve(StopWatch stopWatch) {
        // Passo 1: Converter para a nossa notação: Aqui faz o O ficar do jeito que usamos
        List<List<List<Integer>>> O = new ArrayList<>();

        for (Map<Integer, Integer> map : orders) {
            List<List<Integer>> listOfPairs = new ArrayList<>();
            for (Map.Entry<Integer, Integer> entry : map.entrySet()) {
                listOfPairs.add(Arrays.asList(entry.getKey(), entry.getValue()));
            }
            O.add(listOfPairs);
        }

        // Passo 2: Determinar o NO (número de pedidos)
        int NO = O.size();
        int NI = nItems;





        // Passo 3: Determinar o Io
        List<List<Integer>> Io = buildIo(O, nItems);

        // Imprimir o Io para verificar o resultado
        for (int i = 0; i < Io.size(); i++) {
            System.out.println("Io[" + i + "]: " + Io.get(i));
        }

        
         // Passo 4: Determinar A
List<List<List<Integer>>> A = new ArrayList<>();  // Inicializando A

for (Map<Integer, Integer> map : aisles) {
    List<List<Integer>> listOfPairs = new ArrayList<>();
    for (Map.Entry<Integer, Integer> entry : map.entrySet()) {
        // Adicionando o par (item, capacidade) à lista de pares
        listOfPairs.add(Arrays.asList(entry.getKey(), entry.getValue()));
    }
    // Adicionando a lista de pares ao corredor
    A.add(listOfPairs);
}

int NA = A.size();
int LB = waveSizeLB;
int UB = waveSizeUB;

// Passo 5: Inicializando Ia
List<List<Integer>> Ia = new ArrayList<>();  // Inicializando a lista de listas de corredores para cada item

// Inicializa a lista de corredores para cada item
for (int i = 0; i < NI; i++) {
    Ia.add(new ArrayList<>());  // Cada item começa com uma lista vazia de corredores
}

// Preenchendo Ia com os corredores que possuem cada item
for (int a = 0; a < A.size(); a++) {
    List<List<Integer>> aisle = A.get(a);  // Pega o corredor a (cada corredor é uma lista de pares [item, capacidade])

    for (List<Integer> pair : aisle) {
        int item = pair.get(0);  // Pega o item (primeiro valor do par)
        Ia.get(item).add(a);  // Adiciona o corredor a à lista do item correspondente
    }
}

// Imprimir Ia para verificar o resultado
System.out.println("Corredores Ia:");
for (int i = 0; i < Ia.size(); i++) {
    System.out.println("Ia[" + i + "]: " + Ia.get(i));
}

// Passo 6: Gerando o vetor F
int[] F = new int[NO];
for (int i = 0; i < NO; i++) {
    int soma = 0;
    for (List<Integer> pair : O.get(i)) {
        soma += pair.get(1);
    }
    F[i] = soma;
}

//Passo 6:gera d----------------------------------------------------------------------------------------------

     int[] d = new int[NI]; // Vetor para armazenar a demanda total de cada item

        // Para cada item (i)
        for (int i = 0; i < NI; i++) {
            List<Integer> Io_aux = Io.get(i); // Todos os pedidos que requerem o item i
            int soma_d = 0; // Variável para acumular a soma da demanda de cada item

            // Para cada pedido j que requer o item i
            for (int j = 0; j < Io_aux.size(); j++) {
                int pedidoIndex = Io_aux.get(j); // Índice do pedido
                int sd = buscarValor(O.get(pedidoIndex), i); // Valor da demanda para o item i no pedido
                soma_d += sd; // Soma a demanda
            }

            // Atribui a soma da demanda do item i ao vetor d
            d[i] = soma_d;
        }


// Passo 8: Preprocessamento dos corredores--------------------------------------------------------------------
for (int i = 0; i < NI; i++) { // Para todos os itens
    List<Integer> Ia_aux = Ia.get(i); // Todos os corredores que contêm o item i

    for (int j = 0; j < Ia_aux.size(); j++) {
        int corredor = Ia_aux.get(j); // Índice do corredor
        List<List<Integer>> aisle = A.get(corredor); // Obtém a lista de pares do corredor

        for (List<Integer> pair : aisle) {
            if (pair.get(0) == i && pair.get(1) > d[i]) { // Se a capacidade do item no corredor for maior que d[i]
                pair.set(1, d[i]); // Ajusta a capacidade na variável A
            }
        }
    }
}
// Passo 8: Remover corredores repetidos
for (int i = 0; i < NA - 1; i++) { // Para cada corredor
    for (int j = i + 1; j < NA; j++) { // Compara com todos os corredores seguintes
        if (A.get(i).equals(A.get(j))) { // Se o corredor j é igual ao corredor i
            A.set(j, new ArrayList<>(Arrays.asList(Arrays.asList(1, -1)))); // Define como corredor removido (substitui com [1, -1])
        }
    }
}





 // Passo 9: Executar o novo pré-processamento
    PreprocessingResult preprocessResult = preprocessing(
        NO,         // Número de pedidos (já calculado)
        NI,         // Número de itens (nItems)
        NA,         // Número original de corredores (A.size())
        O,          // Lista de pedidos convertida (List<List<List<Integer>>>)
        A,          // Lista de corredores convertida (List<List<List<Integer>>>)
        waveSizeLB, // Lower bound do tamanho da onda (int)
        waveSizeUB, // Upper bound do tamanho da onda (int)
        Io,         // Lista de pedidos por item (List<List<Integer>>)
        Ia,         // Lista de corredores por item (List<List<Integer>>)
        d           // Vetor de demandas (int[])
    );


    // Passo 10: Atualizar estruturas com os resultados
    A = preprocessResult.A;       // Corredores atualizados (com [1, -1] nos inválidos)
    d = preprocessResult.d;       // Demandas (inalteradas)
    int[] CAP = preprocessResult.CAP; // Novas capacidades
    List<Set<Integer>> ITEM_A = preprocessResult.ITEM_A; // Itens válidos por corredor
    List<List<Integer>> G = preprocessResult.G;         // Grupos de pedidos viáveis
    List<List<Integer>> fff = preprocessResult.fff;         // Grupos de pedidos viáveis

//  Preprocessamento dos corredores--------------------------------------------------------------------
for (int i = 0; i < NI; i++) { // Para todos os itens
    List<Integer> Ia_aux = Ia.get(i); // Todos os corredores que contêm o item i

    for (int j = 0; j < Ia_aux.size(); j++) {
        int corredor = Ia_aux.get(j); // Índice do corredor
        List<List<Integer>> aisle = A.get(corredor); // Obtém a lista de pares do corredor

        for (List<Integer> pair : aisle) {
            if (pair.get(0) == i && pair.get(1) > d[i]) { // Se a capacidade do item no corredor for maior que d[i]
                pair.set(1, d[i]); // Ajusta a capacidade na variável A
            }
        }
    }
}


// Passo 11: Calcular CI e processar dados
        List<Integer> CI = menoresIndices(CAP, 0.8); // 1 - 0.2 = 80%
        Set<Integer> auxSet = new HashSet<>(preprocessResult.aux);
        CI.addAll(auxSet);


        
        // Calcular ITEM_N
        List<Integer> ITEM_N = new ArrayList<>();
        for (int i = 0; i < Ia.size(); i++) {
            if (Ia.get(i).size() == 1 && CI.contains(Ia.get(i).get(0))) {
                ITEM_N.add(i);
            }
        }
        
        // Calcular O_N
        Set<Integer> O_N = new HashSet<>();
        for (int item : ITEM_N) {
            O_N.addAll(Io.get(item));
        }

        // Passo 12: Imprimir dados(opcional)

// ADICIONE AQUI (logo no início do método, antes da criação do modelo)
    System.out.println("\n=== DADOS DE ENTRADA ===");
    
    System.out.println("\n=== INÍCIO DA EXECUÇÃO ===");



  // Passo 13: Executar a heuristica
  // Buscar solução para valores crescentes de L (Primeira iteração)
    int L = 1;
    int L_INICIAL = -1;  // Inicializado com -1 (valor inválido para indicar que ainda não foi encontrado)
    double CB = 0.0;  // Alterado para double para manter precisão
    int maxL = NA;
    List<Integer> bestOrders = new ArrayList<>();
    List<Integer> bestAisles = new ArrayList<>();
   // ============= CONFIGURAÇÃO DE TEMPO =============
    final long SAFETY_MARGIN = 5000; // 5 segundos em milissegundos
    final long EFFECTIVE_MAX_TIME = MAX_RUNTIME - SAFETY_MARGIN;

    // ============= PRIMEIRA FASE (L CRESCENTE) =============
    
    while (L <= maxL && L <= NA) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        // Verificação de tempo crítico
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 1] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        // Atualizar L_INICIAL
        if (L_INICIAL == -1 && result.getObjectiveValue() > 0) {
            L_INICIAL = L;
            System.out.println("[FASE 1] Primeira solução factível em L=" + L);
        }

        // Atualizar melhor solução
        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        // Verificação pós-execução
        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 1] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L++;
    }

    // ============= SEGUNDA FASE (L DECRESCENTE) =============
    CI = menoresIndices(CAP, 0.6);
    CI.addAll(new HashSet<>(preprocessResult.aux));
    
    ITEM_N.clear();
    for (int i = 0; i < Ia.size(); i++) {
        if (Ia.get(i).size() == 1 && CI.contains(Ia.get(i).get(0))) {
            ITEM_N.add(i);
        }
    }
    
    O_N.clear();
    for (int item : ITEM_N) {
        O_N.addAll(Io.get(item));
    }

    L = (L_INICIAL > 0) ? L_INICIAL - 1 : 1;
    while (L > 0) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 2] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if (result.getObjectiveValue() > 0) {
            L_INICIAL = Math.min(L_INICIAL, L);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 2] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L--;
    }

    // ============= SEGUNDA FASE (L CRESCENTE OTIMIZADO) =============
    L = L_INICIAL;
    
    while (L <= maxL && L <= NA) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 3] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 3] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L++;
    }




// ============= TERCEIRA FASE (L DECRESCENTE) =============
    CI = menoresIndices(CAP, 0.4);
    CI.addAll(new HashSet<>(preprocessResult.aux));
    
    ITEM_N.clear();
    for (int i = 0; i < Ia.size(); i++) {
        if (Ia.get(i).size() == 1 && CI.contains(Ia.get(i).get(0))) {
            ITEM_N.add(i);
        }
    }
    
    O_N.clear();
    for (int item : ITEM_N) {
        O_N.addAll(Io.get(item));
    }

    L = (L_INICIAL > 0) ? L_INICIAL - 1 : 1;
    
    while (L > 0) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 2] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if (result.getObjectiveValue() > 0) {
            L_INICIAL = Math.min(L_INICIAL, L);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 2] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L--;
    }

    // ============= TERCEIRA FASE (L CRESCENTE OTIMIZADO) =============
    L = L_INICIAL;
    
    while (L <= maxL && L <= NA) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 3] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 3] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L++;
    }





// ============= QUARTA FASE (L DECRESCENTE) =============
    CI = menoresIndices(CAP, 0.2);
    CI.addAll(new HashSet<>(preprocessResult.aux));
    
    ITEM_N.clear();
    for (int i = 0; i < Ia.size(); i++) {
        if (Ia.get(i).size() == 1 && CI.contains(Ia.get(i).get(0))) {
            ITEM_N.add(i);
        }
    }
    
    O_N.clear();
    for (int item : ITEM_N) {
        O_N.addAll(Io.get(item));
    }

    L = (L_INICIAL > 0) ? L_INICIAL - 1 : 1;
    
    while (L > 0) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 2] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if (result.getObjectiveValue() > 0) {
            L_INICIAL = Math.min(L_INICIAL, L);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 2] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L--;
    }

    // ============= QUARTA FASE (L CRESCENTE OTIMIZADO) =============
    L = (L_INICIAL > 0) ? L_INICIAL - 1 : 1;
    
    while (L <= maxL && L <= NA) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 3] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 3] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L++;
    }



// ============= ULTIMA FASE (L DECRESCENTE) =============
CI.clear();  // Limpa a lista existente
CI.addAll(preprocessResult.aux);  // Adiciona todos elementos do Set
        


    ITEM_N.clear();
    for (int i = 0; i < Ia.size(); i++) {
        if (Ia.get(i).size() == 1 && CI.contains(Ia.get(i).get(0))) {
            ITEM_N.add(i);
        }
    }
    
    O_N.clear();
    for (int item : ITEM_N) {
        O_N.addAll(Io.get(item));
    }

    L = (L_INICIAL > 0) ? L_INICIAL - 1 : 1;
    
    while (L > 0) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 2] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if (result.getObjectiveValue() > 0) {
            L_INICIAL = Math.min(L_INICIAL, L);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 2] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L--;
    }

    // ============= TERCEIRA FASE (L CRESCENTE OTIMIZADO) =============
    L = L_INICIAL;
    
    while (L <= maxL && L <= NA) {
        long remainingTimeMs = EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS);
        
        if (remainingTimeMs < 1000) {
            System.out.println("[FASE 3] Tempo crítico - retornando solução atual");
            break;
        }
        
        OptimizationResult result = runModel(
            (int) Math.ceil(CB), LB, UB, NO, NA, NI, L,
            F, A, Ia, O, Io, CI, G, fff,
            ITEM_N, O_N, remainingTimeMs / 1000.0
        );

        if (result.getObjectiveValue() > CB) {
            CB = result.getObjectiveValue();
            bestOrders = new ArrayList<>(result.getSelectedOrders());
            bestAisles = new ArrayList<>(result.getSelectedAisles());
            maxL = Math.min((int) (UB / CB), NA);
        }

        if ((EFFECTIVE_MAX_TIME - stopWatch.getTime(TimeUnit.MILLISECONDS)) < 1000) {
            System.out.println("[FASE 3] Tempo crítico pós-execucao - saindo");
            break;
        }
        
        L++;
    }







    System.out.println("Tempo total usado: " + stopWatch.getTime(TimeUnit.SECONDS) + "s");
    return new ChallengeSolution(
        new HashSet<>(bestOrders),
        new HashSet<>(bestAisles)
    );




 // Acho que está ok, conferir dpz se precisa atualizar Ia, neste caso tomar cuidado com Ia duplicado.
//
//
//
    // Passo 7: Recalcular Ia
   // List<List<Integer>> newIa = new ArrayList<>();
   // for (int i = 0; i < NI; i++) newIa.add(new ArrayList<>());
   // for (int a = 0; a < A.size(); a++) {
       // if (isInvalidAisle(A.get(a))) continue;
       // for (List<Integer> pair : A.get(a)) {
     //       int item = pair.get(0);
   //         if (item < NI) newIa.get(item).add(a);
  //     }
 //   }
//    Ia = newIa;

    }




public class HeuristicSolver {
    
    public static class HeuristicResult {
        public double bestObjective;  // Melhor valor FO encontrado
        public List<Integer> selectedAisles;  // Corredores selecionados (A_solution)
        public List<Integer> selectedOrders;   // Pedidos selecionados (O_solution)
        
        public HeuristicResult(double bestObjective, 
                             List<Integer> selectedAisles, 
                             List<Integer> selectedOrders) {
            this.bestObjective = bestObjective;
            this.selectedAisles = selectedAisles;
            this.selectedOrders = selectedOrders;
        }
    }

    public HeuristicResult runHeuristic(double initialFO, double timeLimit, int initialL, 
                                      int numOrders, int numItems, int numAisles,
                                      List<List<List<Integer>>> orders, 
                                      List<List<List<Integer>>> aisles,
                                      int lowerBound, int upperBound,
                                      List<List<Integer>> itemToOrdersMap,
                                      List<List<Integer>> itemToAislesMap,
                                      List<Integer> fixedAisles,
                                      List<List<Integer>> orderGroups,
                                      List<List<Integer>> additionalConstraints,
                                      List<Integer> restrictedItems,
                                      Set<Integer> excludedOrders,
                                      StopWatch timer) {
        
        double bestFO = initialFO;
        List<Integer> bestAisles = new ArrayList<>();
        List<Integer> bestOrders = new ArrayList<>();
        double maxL = upperBound / initialFO;
        int currentL = initialL;

        // Primeira fase: testa valores de L decrescentes
        while (currentL > 0) {
            if (timer.getTime(TimeUnit.SECONDS) >= timeLimit) {
                System.out.println("Tempo esgotado - Fase decrescente");
                break;
            }

            OptimizationResult result = solveModel((int)initialFO, lowerBound, upperBound, 
                                                 numOrders, numAisles, numItems, currentL,
                                                 orders, aisles, itemToOrdersMap, itemToAislesMap,
                                                 fixedAisles, orderGroups, additionalConstraints,
                                                 restrictedItems, excludedOrders, 
                                                 timeLimit - timer.getTime(TimeUnit.SECONDS));

            // Atualiza a melhor solução se encontrou uma melhor
            if (result.getObjectiveValue() > bestFO) {
                bestFO = result.getObjectiveValue();
                bestAisles = result.getSelectedAisles();
                bestOrders = result.getSelectedOrders();
                maxL = upperBound / bestFO;
            }

            currentL--;
            if (result.getObjectiveValue() == 0) {
                break;
            }
        }

        // Segunda fase: testa valores de L crescentes
        currentL = initialL + 1;
        while (currentL <= maxL && currentL <= numAisles) {
            if (timer.getTime(TimeUnit.SECONDS) >= timeLimit) {
                System.out.println("Tempo esgotado - Fase crescente");
                break;
            }

            OptimizationResult result = solveModel((int)initialFO, lowerBound, upperBound,
                                                numOrders, numAisles, numItems, currentL,
                                                orders, aisles, itemToOrdersMap, itemToAislesMap,
                                                fixedAisles, orderGroups, additionalConstraints,
                                                restrictedItems, excludedOrders,
                                                timeLimit - timer.getTime(TimeUnit.SECONDS));

            // Atualiza a melhor solução se encontrou uma melhor
            if (result.getObjectiveValue() > bestFO) {
                bestFO = result.getObjectiveValue();
                bestAisles = result.getSelectedAisles();
                bestOrders = result.getSelectedOrders();
                maxL = upperBound / bestFO;
            }

            currentL++;
        }

        return new HeuristicResult(bestFO, bestAisles, bestOrders);
    }

    private OptimizationResult solveModel(int currentBest, int lb, int ub,
                                       int no, int na, int ni, int l,
                                       List<List<List<Integer>>> O,
                                       List<List<List<Integer>>> A,
                                       List<List<Integer>> Io,
                                       List<List<Integer>> Ia,
                                       List<Integer> CI,
                                       List<List<Integer>> G,
                                       List<List<Integer>> fff,
                                       List<Integer> ITEM_N,
                                       Set<Integer> O_N,
                                       double remainingTime) {
        // Implementação do seu solver de otimização
        // Retorna um OptimizationResult contendo:
        // - Valor objetivo
        // - Lista de corredores selecionados
        // - Lista de pedidos selecionados
        return new OptimizationResult(0.0, new ArrayList<>(), new ArrayList<>());
    }
}









public static List<Integer> menoresIndices(int[] maxC, double frac) {
    List<Integer> indices = new ArrayList<>();
    for (int i = 0; i < maxC.length; i++) indices.add(i);
    Collections.sort(indices, Comparator.comparingInt(i -> maxC[i]));
    int k = Math.max(0, Math.min((int) Math.floor(frac * maxC.length), maxC.length));
    return indices.subList(0, k); 
}















private void printPreprocessingData(PreprocessingResult result) {
    System.out.println("\n========== DADOS DO PRÉ-PROCESSAMENTO ==========");
    
    // 1. Imprimir A (corredores)
    System.out.println("\nA (Corredores):");
    for (int i = 0; i < result.A.size(); i++) {
        System.out.println("Corredor " + i + ": " + result.A.get(i));
    }
    
    // 2. Imprimir d (demandas)
    System.out.println("\nd (Demandas): " + Arrays.toString(result.d));
    
    // 3. Imprimir CAP (capacidades)
    System.out.println("\nCAP (Capacidades): " + Arrays.toString(result.CAP));
    
    // 4. Imprimir ITEM_A (itens com capacidade suficiente)
    System.out.println("\nITEM_A:");
    for (int i = 0; i < result.ITEM_A.size(); i++) {
        System.out.println("Corredor " + i + ": " + result.ITEM_A.get(i));
    }
    
    // 5. Imprimir ITEM_B (todos os itens)
    System.out.println("\nITEM_B:");
    for (int i = 0; i < result.ITEM_B.size(); i++) {
        System.out.println("Corredor " + i + ": " + result.ITEM_B.get(i));
    }
    
    // 6. Imprimir DOM (subconjuntos)
    System.out.println("\nDOM (Subconjuntos):");
    for (int i = 0; i < result.DOM.size(); i++) {
        System.out.println("Corredor " + i + " domina: " + result.DOM.get(i));
    }
    
    // 7. Imprimir aux (corredores marcados)
    System.out.println("\naux (Corredores Invalidos): " + result.aux);
    
    // 8. Imprimir ITEM_O (itens por pedido)
    System.out.println("\nITEM_O (Itens por Pedido):");
    for (int i = 0; i < result.ITEM_O.size(); i++) {
        System.out.println("Pedido " + i + ": " + result.ITEM_O.get(i));
    }
    
    // 9. Imprimir fff (pedidos que cabem nos corredores)
    System.out.println("\nfff:");
    for (int i = 0; i < result.fff.size(); i++) {
        System.out.println("Corredor " + i + ": " + result.fff.get(i));
    }
    
    // 10. Imprimir ggg (pedidos exclusivos)
    System.out.println("\nggg:");
    for (int i = 0; i < result.ggg.size(); i++) {
        System.out.println("Corredor " + i + ": " + result.ggg.get(i));
    }
    
    // 11. Imprimir G (grupos viáveis)
    System.out.println("\nG (Grupos Viáveis):");
    for (int i = 0; i < result.G.size(); i++) {
        System.out.println("Corredor " + i + ": " + result.G.get(i));
    }
    
    System.out.println("==============================================\n");
}


    // Método para buscar o valor de um item em uma lista de pedidos
    private int buscarValor(List<List<Integer>> pedidos, int item) {
        int valor = 0;
        for (List<Integer> pedido : pedidos) {
            int pedidoItem = pedido.get(0);  
            int quantidade = pedido.get(1);  

            if (pedidoItem == item) {
                valor = quantidade;  
                break;  
            }
        }
        return valor;  
    }

private PreprocessingResult preprocessing(
    int NO, int NI, int NA, 
    List<List<List<Integer>>> O, 
    List<List<List<Integer>>> A,
    int LB, int UB,
    List<List<Integer>> Io,
    List<List<Integer>> Ia,
    int[] d) {

    // Passo 1: Calcular ITEM_A e ITEM_B
    List<Set<Integer>> ITEM_A = new ArrayList<>();
    List<Set<Integer>> ITEM_B = new ArrayList<>();
    for (List<List<Integer>> aisle : A) {
        Set<Integer> aItems = new HashSet<>();
        Set<Integer> bItems = new HashSet<>();
        for (List<Integer> pair : aisle) {
            int item = pair.get(0);
            bItems.add(item);
            if (pair.get(1) >= d[item]) {
                aItems.add(item);
            }
        }
        ITEM_A.add(aItems);
        ITEM_B.add(bItems);
    }

    // Passo 2: Calcular DOM
    List<List<Integer>> DOM = new ArrayList<>();
    for (int j = 0; j < ITEM_A.size(); j++) {
        List<Integer> indices = new ArrayList<>();
        Set<Integer> baseSet = ITEM_A.get(j);
        for (int i = 0; i < ITEM_B.size(); i++) {
            if (i != j && baseSet.containsAll(ITEM_B.get(i))) {
                indices.add(i);
            }
        }
        DOM.add(indices);
    }

    // Passo 3: Calcular aux
    Set<Integer> aux = new HashSet<>();
    for (List<Integer> list : DOM) aux.addAll(list);
    
    Set<Integer> aux2 = new HashSet<>();
    for (int i = 0; i < A.size(); i++) {
        for (List<Integer> pair : A.get(i)) {
            if (pair.get(0) == 1 && pair.get(1) == -1) {
                aux2.add(i);
                break;
            }
        }
    }
    aux.addAll(aux2);

    // Passo 4: Marcar corredores inválidos
for (Integer invalidIdx : aux) {
    if (invalidIdx < A.size()) {
        A.set(invalidIdx, new ArrayList<>(Arrays.asList(Arrays.asList(1, -1))));
        //                                                              ↑↑ Correção aqui
    }
}

    // Passo 5: Calcular ITEM_O
    List<Set<Integer>> ITEM_O = new ArrayList<>();
    for (List<List<Integer>> order : O) {
        ITEM_O.add(order.stream()
            .map(pair -> pair.get(0))
            .collect(Collectors.toSet()));
    }

    // Passo 6: Calcular fff
    List<List<Integer>> fff = new ArrayList<>();
    for (int i = 0; i < A.size(); i++) {
        List<Integer> indices = new ArrayList<>();
        if (!isInvalidAisle(A.get(i))) {
            Set<Integer> aItems = ITEM_A.get(i);
            for (int j = 0; j < ITEM_O.size(); j++) {
                if (aItems.containsAll(ITEM_O.get(j))) {
                    indices.add(j);
                }
            }
        }
        fff.add(indices);
    }

    // Passo 7: Recalcular Ia
    List<List<Integer>> newIa = new ArrayList<>();
    for (int i = 0; i < NI; i++) newIa.add(new ArrayList<>());
    for (int a = 0; a < A.size(); a++) {
        if (isInvalidAisle(A.get(a))) continue;
        for (List<Integer> pair : A.get(a)) {
            int item = pair.get(0);
            if (item < NI) newIa.get(item).add(a);
        }
    }
    Ia = newIa;

    // Passo 8: Calcular ggg
    List<List<Integer>> ggg = new ArrayList<>();
    for (int i = 0; i < NA; i++) ggg.add(new ArrayList<>());
    for (int item = 0; item < NI; item++) {
        if (Ia.get(item).size() == 1) {
            int aisleIdx = Ia.get(item).get(0);
            if (aisleIdx < NA && !Io.get(item).isEmpty()) {
                ggg.get(aisleIdx).add(Io.get(item).get(0));
            }
        }
    }

    // Passo 9: Calcular G
    List<List<Integer>> G = new ArrayList<>();
    for (int i = 0; i < NA; i++) {
        Set<Integer> intersection = new HashSet<>(fff.get(i));
        intersection.retainAll(new HashSet<>(ggg.get(i)));
        G.add(new ArrayList<>(intersection));
    }


// Passo 10: Recalcular CAP garantindo CAP <= d
int[] CAP = new int[NA];
for (int a = 0; a < NA; a++) {
    if (isInvalidAisle(A.get(a))) {
        CAP[a] = -1;
    } else {
        // Calcula a capacidade original do corredor
        int originalCap = A.get(a).stream()
                          .mapToInt(p -> p.get(1))
                          .sum();
        
        // Para cada item no corredor, verifica a demanda
        int limitedCap = originalCap;
        for (List<Integer> item : A.get(a)) {
            int itemId = item.get(0);
            int itemCap = item.get(1);
            
            // Se a capacidade do item no corredor > demanda total, ajusta
            if (itemCap > d[itemId]) {
                limitedCap = limitedCap - itemCap + d[itemId];
            }
        }
        
        CAP[a] = limitedCap;
    }
}
    return new PreprocessingResult(
        A, d, CAP, ITEM_A, ITEM_B, DOM, aux, ITEM_O, fff, ggg, G
    );
}

private boolean isInvalidAisle(List<List<Integer>> aisle) {
    return aisle.size() == 1 && 
           aisle.get(0).get(0) == 1 && 
           aisle.get(0).get(1) == -1;
}

    // Método para construir a lista Io
    public static List<List<Integer>> buildIo(List<List<List<Integer>>> O, int nItems) {
        // Inicializa Io com nItems listas vazias
        List<List<Integer>> Io = new ArrayList<>();
        for (int i = 0; i < nItems; i++) {
            Io.add(new ArrayList<>());  // Cada tipo de item recebe uma lista vazia
        }

        // Preenche Io conforme os pedidos em O
        for (int o = 0; o < O.size(); o++) {
            List<List<Integer>> order = O.get(o);  // Pega o pedido o
            for (List<Integer> itemPair : order) {
                int item = itemPair.get(0);  // Pega o item (primeiro valor do par)
                Io.get(item).add(o);  // Adiciona o pedido o à lista correspondente ao item
            }
        }

        return Io;  // Retorna a lista Io preenchida
    }

    /*
     * Get the remaining time in seconds
     */
    protected long getRemainingTime(StopWatch stopWatch) {
        return Math.max(
                TimeUnit.SECONDS.convert(MAX_RUNTIME - stopWatch.getTime(TimeUnit.MILLISECONDS), TimeUnit.MILLISECONDS),
                0);
    }

    protected boolean isSolutionFeasible(ChallengeSolution challengeSolution) {
        Set<Integer> selectedOrders = challengeSolution.orders();
        Set<Integer> visitedAisles = challengeSolution.aisles();
        if (selectedOrders == null || visitedAisles == null || selectedOrders.isEmpty() || visitedAisles.isEmpty()) {
            return false;
        }

        int[] totalUnitsPicked = new int[nItems];
        int[] totalUnitsAvailable = new int[nItems];

        // Calculate total units picked
        for (int order : selectedOrders) {
            for (Map.Entry<Integer, Integer> entry : orders.get(order).entrySet()) {
                totalUnitsPicked[entry.getKey()] += entry.getValue();
            }
        }

        // Calculate total units available
        for (int aisle : visitedAisles) {
            for (Map.Entry<Integer, Integer> entry : aisles.get(aisle).entrySet()) {
                totalUnitsAvailable[entry.getKey()] += entry.getValue();
            }
        }

        // Check if the total units picked are within bounds
        int totalUnits = Arrays.stream(totalUnitsPicked).sum();
        if (totalUnits < waveSizeLB || totalUnits > waveSizeUB) {
            return false;
        }

        // Check if the units picked do not exceed the units available
        for (int i = 0; i < nItems; i++) {
            if (totalUnitsPicked[i] > totalUnitsAvailable[i]) {
                return false;
            }
        }

        return true;
    }

    protected double computeObjectiveFunction(ChallengeSolution challengeSolution) {
        Set<Integer> selectedOrders = challengeSolution.orders();
        Set<Integer> visitedAisles = challengeSolution.aisles();
        if (selectedOrders == null || visitedAisles == null || selectedOrders.isEmpty() || visitedAisles.isEmpty()) {
            return 0.0;
        }
        int totalUnitsPicked = 0;

        // Calculate total units picked
        for (int order : selectedOrders) {
            totalUnitsPicked += orders.get(order).values().stream()
                    .mapToInt(Integer::intValue)
                    .sum();
        }

        // Calculate the number of visited aisles
        int numVisitedAisles = visitedAisles.size();

        // Objective function: total units picked / number of visited aisles
        return (double) totalUnitsPicked / numVisitedAisles;
    }
}
