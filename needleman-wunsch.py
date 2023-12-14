import numpy as np
import math

# Pontuações configuradas conforme pedido
MATCH_SCORE = +7 # Ganho de pontuação por paridade de proteínas
MISMATCH_SCORE = -3 # Perda de pontuação por disparidade de proteínas
GAP_SCORE = -4 # Perda de pontuação por burado no pareamento

# Primeira sequência para o alinhamento
SEQUENCE_ONE = 'GCATT'
# Segunda sequência para o alinhamento
SEQUENCE_TWO = 'GTT'

# Primeira sequência alinhada
alignment_a = ''
# Segunda sequência alinhada
alignment_b = ''

# Enumeradores para escolha de caminho na matrix
DIG = 1
LEFT = 2
UP = 3

# Quantidade de colunas baseadas na primeira sequência 
columns = len(SEQUENCE_ONE) + 1
# Quantidade de linhas baseadas na segunda sequência
rows = len(SEQUENCE_TWO) + 1

# Matriz de armazenamento de pontuação
matrix = np.zeros((columns, rows))
# Matriz de armazenemento de caminho
matrix_path = np.zeros((columns, rows))

# Preenchimento de bordas baseado em preenchimento
for i in range (1, columns):
    matrix[i][0] = i * GAP_SCORE

# Preenchimento de bordas baseado em preenchimento
for j in range (1, rows):
    matrix[0][j] = j * GAP_SCORE

# População da matriz com pontuação e caminho
def matrix_population(i, j):
    # Percorre inversamento, populando os nós dependentes primeiro
    # # O ponto (5,5) da matriz depende dos pontos (5,4), (4,5) e (4,4) para seu calculo
    # # Logo a recursão popula primeiro os pontos (5,4), (4,5) e (4,4)
    if i > 1 and matrix[i - 1][j] == 0: # Esquerda
        matrix_population(i - 1, j)
    if i > 1 and j > 1 and matrix[i - 1][j - 1] == 0: # Diagonal
        matrix_population(i - 1, j - 1)
    if j > 1 and matrix[i][j - 1] == 0: # Cima
        matrix_population(i, j - 1)

    # Letra da primeira sequência
    letter_one = SEQUENCE_ONE[i - 1]
    # Letra da segunda sequência
    letter_two = SEQUENCE_TWO[j - 1]

    # Pontuação do ponto da diagonal
    dig = matrix[i - 1][j - 1]
    # Pontuação do ponto da esquerda
    left = matrix[i][j - 1]
    # Pontuação do ponto de cima
    up = matrix[i - 1][j]

    # Pontuação do ponto atual considerando um deslocamento da diagonal
    score_from_dig = dig + (MATCH_SCORE if letter_one == letter_two else MISMATCH_SCORE)
    # Pontuação do ponto atual considerando um deslocamento da esquerda
    score_from_left = left + GAP_SCORE
    # Pontuação do ponto atual considerando um deslocamento de cima
    score_from_up = up + GAP_SCORE

    # Escolhe a melhor opção entre diagonal, esquerda e cima
    new_score = max(score_from_dig, score_from_up, score_from_left)

    # Adiciona a pontuação ao ponto atual
    matrix[i][j] = new_score

    # Múltiplos deslocamentos podem possuir a mesma pontuação
    # Escolhido o caminho com a menor pontuação anterior

    # Inicializa os valores em menos infinito
    score_dig = -math.inf
    score_left = -math.inf
    score_up = -math.inf

    # Caso a pontuação selecionada tenha vindo da diagonal
    if new_score == score_from_dig:
        # Salva a pontuação da diagonal na variável score_dig
        score_dig = matrix[i - 1][j - 1]
        
    # Caso a pontuação selecionada tenha vindo de cima
    if new_score == score_from_up:
        # Salva a pontuação de cima na variável score_up
        score_up = matrix[i - 1][j]
    
    # Caso a pontuação selecionada tenha vindo da esquerda
    if new_score == score_from_left:
        # Salva a pontuação da esquerda na variável score_left
        score_left = matrix[i][j - 1]

    # Escolhe a melhor pontuação anterior
    best_score = max(score_dig, score_up, score_left)

    # Salva o deslocamento na matriz matrix_path
    if best_score == score_dig:
        matrix_path[i][j] = DIG
    elif best_score == score_up:
        matrix_path[i][j] = UP
    else:
        matrix_path[i][j] = LEFT

# Leitura do caminho com criação dos alinhamento
def matrix_read(i, j):
    global alignment_a
    global alignment_b

    # Para ao atingir o primeiro ponto da matriz
    if i == 0 and j == 0:
        return

    # Lê a escolha baseado na matriz
    choice = matrix_path[i][j]
    
    if choice == DIG: # Se a escolha for diagonal
        # O resultado é paridade, logo adiciona as letras das sequências (iguais) aos alinhamentos
        alignment_a = SEQUENCE_ONE[i - 1] + alignment_a
        alignment_b = SEQUENCE_TWO[j - 1] + alignment_b
        matrix_read(i - 1, j - 1) # Continua a recursão no deslocamento escolhido
    elif choice == UP: # Se a escolha for cima
        # O resultado foi um preenchimento na segunda sequência
        alignment_a = SEQUENCE_ONE[i - 1] + alignment_a # Adiciona a letra na primeira sequência
        alignment_b = '-' + alignment_b # Adiciona um preenchimento na segunda sequência
        matrix_read(i - 1, j) # Continua a recursão no deslocamento escolhido
    elif choice == LEFT: # Se a escolha for esquerda
        # O resultado foi um preenchimento na primeira sequência
        alignment_a = '-' + alignment_a  # Adiciona um preenchimento na primeira sequência
        alignment_b = SEQUENCE_TWO[j - 1] + alignment_b # Adiciona a letra na segunda sequência
        matrix_read(i, j - 1) # Continua a recursão no deslocamento escolhido

# Inicia a população dos dados do ponto mais inferior da direita
matrix_population(columns - 1, rows - 1)

# Inicia a leitura dos alinhamentos do ponto mais inferior da direita
matrix_read(columns - 1, rows - 1)

# Impressão de informações no console
print('Score: ' + str(matrix[columns - 1][rows - 1]))
print(alignment_a)
print(alignment_b)

print(matrix)
print(matrix_path)