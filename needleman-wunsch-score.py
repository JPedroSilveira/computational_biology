import time
import sys

# Pontuações configuradas conforme pedido
MATCH_SCORE = +1 # Ganho de pontuação por paridade de proteínas
MISMATCH_SCORE = -1 # Perda de pontuação por disparidade de proteínas
GAP_SCORE = -2 # Perda de pontuação por burado no pareamento

# Primeira sequência para o alinhamento
SEQUENCE_ONE = ''
# Segunda sequência para o alinhamento
SEQUENCE_TWO = ''

# Número das sequencias deve ser passados logo após a chamada do script como argumentos
# Examplo para comprar sequências 6 e 7: `python3 needleman-wunsch-score.py 6 7`
FIRST_SEQUENCE_NAME = sys.argv[1]
SECOND_SEQUENCE_NAME = sys.argv[2]

# Define um limite de tamanho para analisar cada sequência
# Motivação: limitação no tempo de execução e memória
# Exemplo: dada uma sequência de 80.000 caracteres e SEQUENCE_LIMIT = 40000, 
#   serão análisados apenas do primeiro caractere até o que estiver na posição 40000
SEQUENCE_LIMIT = 40000

# Leitura da primeira sequência
with open('dataset_seq_two/seq_' + FIRST_SEQUENCE_NAME + '.fasta.fna', 'r') as file:
    SEQUENCE_ONE = file.read().replace('\n', '')

# Leitura da segunda sequência
with open('dataset_seq_two/seq_' + SECOND_SEQUENCE_NAME + '.fasta.fna', 'r') as file:
    SEQUENCE_TWO = file.read().replace('\n', '')

# Informações de execução
print('Executing with sequences ' + FIRST_SEQUENCE_NAME + ' and ' + SECOND_SEQUENCE_NAME)
print(' => with sizes ' + str(len(SEQUENCE_ONE)) + ' and ' + str(len(SEQUENCE_TWO))) 
print(' => using a limit of ' + str(SEQUENCE_LIMIT) + '\n')

# Define a sequência conforme a limitação definidas
SEQUENCE_ONE = SEQUENCE_ONE[0:SEQUENCE_LIMIT]
SEQUENCE_TWO = SEQUENCE_TWO[0:SEQUENCE_LIMIT]

# Calcula a pontuação de semelhança entre as sequências
def calculate_score(sequenceOne, sequenceTwo):
    # Define a maior e a menor sequência
    if len(sequenceOne) > len(sequenceTwo):     
        biggerSequence = sequenceOne
        smallerSequence = sequenceTwo
    else:
        biggerSequence = sequenceTwo
        smallerSequence = sequenceOne  

    # Utilitário para armazenas o tamanho de cada sequência somado de 1
    biggerSequenceSizePlusOne = len(biggerSequence) + 1
    smallerSequenceSizePlusOne = len(smallerSequence) + 1  

    # Realiza o cálculo pelas diagonais mantendo em memória apenas os dois últimos vetores em vez da matriz completa
    line_one = [0] # Por definição a primeira diagonal apenas conterá o valor zero, dado que compara a primeira sequência com ela mesma
    line_two = [GAP_SCORE, GAP_SCORE] # A segunda diagonal por definição terá apenas os valores de um gap entre a primeira e a segunda sequência
    diagSize = 1 # Mantem o tamanho da diagonal

    # O tamanho da diagonal aumenta até alcançar a diagonal que passa pela última proteína da menor sequência 
    # O tamanho da diagonal começa a decair após passar pela diagonal que passa última proteína
    # Contador que aumenta para cada iteração do laço até alcançar o tamanho da maior sequência
    countUntilReachBiggerSequenceSize = 0 
    # Contador que aumenta para cada iteração do laço após a quantidade de iterações passar o tamanho da maior sequência 
    countSinceReachBiggerSequenceSizeUntilZero = 0

    # Inicia na terceira diagonal, até a última diagonal possível
    for l in range(3, biggerSequenceSizePlusOne + smallerSequenceSizePlusOne):
        # Inicia as distâncias da próxima diagonal
        new_line = []

        # Preenchimento de borda se necessário
        if biggerSequenceSizePlusOne >= l:
            new_line.append((l - 1) * GAP_SCORE)

        # De 0 até o tamanho da diagonal menos 2 por causa das bordas
        for i in range(0, diagSize):
            # Seleção da pontuação nova
            new_score = calculate_new_score(biggerSequence, smallerSequence, line_one, line_two, i, countUntilReachBiggerSequenceSize, countSinceReachBiggerSequenceSizeUntilZero)

            # Adição a nova linha
            new_line.append(new_score)

        if smallerSequenceSizePlusOne >= l:
            # Preenchimento de borda
            new_line.append((l - 1) * GAP_SCORE)
            # Aumenta o tamanho da diagonal
            diagSize += 1
        elif biggerSequenceSizePlusOne < l:
            # Diminui o tamanho da diagonal
            diagSize -= 1

        # Atualiza o contator baseado no momento de execução
        if countUntilReachBiggerSequenceSize + 1 < len(biggerSequence):
            countUntilReachBiggerSequenceSize += 1
        else:
            countSinceReachBiggerSequenceSizeUntilZero += 1

        # Atualiza as diagonais
        line_one = line_two
        line_two = new_line
    
    return new_line[0]

def calculate_new_score(biggerSequence, smallerSequence, line_one, line_two, i, countUntilReachBiggerSequenceSize, countSinceReachBiggerSequenceSizeUntilZero):
    #   Letra da primeira sequência
    letter_one = biggerSequence[countUntilReachBiggerSequenceSize - i]
    # Letra da segunda sequência
    letter_two = smallerSequence[i + countSinceReachBiggerSequenceSizeUntilZero]

    indexAdjust = 0 if countSinceReachBiggerSequenceSizeUntilZero == 0 else 1

    # Pontuação da diagonal
    dig = line_one[i + indexAdjust]
    # Pontuação esquerda
    left = line_two[i + indexAdjust]
    # Pontuação direita
    right = line_two[i + 1 - indexAdjust]

    # Calculo das possíveis pontuações
    score_from_dig = dig + (MATCH_SCORE if letter_one == letter_two else MISMATCH_SCORE)
    score_from_left = left + GAP_SCORE
    score_from_right = right + GAP_SCORE

    # Seleção da pontuação maior
    return max(score_from_dig, score_from_left, score_from_right)

start_time = time.time()

# Inicia a busca pela pontuação
score = calculate_score(SEQUENCE_ONE, SEQUENCE_TWO)

print("Execution time: %s seconds \n" % (time.time() - start_time))

# Impressão do resultado
print('Score: ' + str(score))