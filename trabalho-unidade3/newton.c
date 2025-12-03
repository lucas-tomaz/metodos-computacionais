/*
UNIVERSIDADE FEDERAL DO RIO GRANDE DO NORTE
CENTRO DE TECNOLOGIA - DEPARTAMENTO DE ENGENHARIA ELÉTRICA
DCA0304 - MÉTODOS COMPUTACIONAIS EM ENGENHARIA
DOCENTE: PAULO SERGIO DA MOTTA PIRES
DISCENTES: LUCAS TOMAZ DE MOURA, SOFIA SEVERO GALVÃO, RENATO TAVEIRA DA SILVA E
GUSTAVO HENRIQUE AZEVEDO.

DESCRIÇÃO DO PROGRAMA:
Algoritmo computacional empregado na resolução de sistemas não-lineares
utlilizando o método de Newton-Raphson.

VERSÃO: 1.1.5


*/

#include "math.h"
#include <stdio.h>
#include <time.h>

#define e 2.718281828459045
#define TOL 1e-6
#define MAX_ITER 100

void EliminacaoGauss(float mat[3][4], float a, float b, float c) {
  int i, j, k;
  float pivo;

  // Aplicando a eliminação de Gauss
  for (i = 0; i < 3; i++) {
    // Tornar a diagonal principal igual a 1
    for (j = i + 1; j < 3; j++) {
      pivo = mat[j][i] / mat[i][i];
      for (k = 0; k <= 3; k++) {
        mat[j][k] -= pivo * mat[i][k];
      }
    }
  }

  // Resolvendo o sistema de equações
  printf("Soluções:\n");
  float x[3];
  for (i = 2; i >= 0; i--) {
    x[i] = mat[i][3];
    for (j = i + 1; j < 3; j++) {
      x[i] -= mat[i][j] * x[j];
    }
    x[i] /= mat[i][i];
  }
    float s_1 = x[0] + a;
    float s_2 = x[1] + b;
    float s_3 = x[2] + c;

    printf("s_1 = %.6f\n", s_1);
    printf("s_2 = %.6f\n", s_2);
    printf("s_3 = %.6f\n", s_3);
    
    float f1, f2, f3;

    f1 = s_1 * s_3 - s_3 * pow(e, pow(s_1, 2)) + pow(10, -4);

    f2 = s_1 * (pow(s_1, 2) + pow(s_2, 2)) + pow(s_2, 2) * (s_3 - s_2);
      
    f3 = (pow(s_1, 3) + pow(s_3, 3));
    
    printf("Solução de f_1 = %.6f\n", f1);
    printf("Solução de f_2 = %.6f\n", f2);
    printf("Solução de f_3 = %.6f\n", f3);
}

int main(void) {
  int i, j, iteracoes = 0;
  clock_t inicio, fim;
  double tempo_processamento;

  // VARIÁVEIS DE CONTROLE DO SISTEMA NÃ0-LINEAR(CHUTE INICIAL)
  float x_1 = 0.01;
  float x_2 = 0.01;
  float x_3 = -0.01;

  float f_1, f_2, f_3;
  float norma_f;

  // Iniciar medição de tempo
  inicio = clock();

  // SISTEMA NÃO-LINEAR:
  f_1 = x_1 * x_3 - x_3 * pow(e, pow(x_1, 2)) + pow(10, -4);

  f_2 = x_1 * (pow(x_1, 2) + pow(x_2, 2)) + pow(x_2, 2) * (x_3 - x_2);

  f_3 = (pow(x_1, 3) + pow(x_3, 3));

  // JACOBIANA DO SISTEMA NÃO-LINEAR:
  float jacobiana[3][3];
  jacobiana[0][0] = x_3 - 2 * x_1 * x_3 * pow(e, pow(x_1, 2));
  jacobiana[0][1] = 0;
  jacobiana[0][2] = x_1 - pow(e, pow(x_1, 2));

  jacobiana[1][0] = 3 * pow(x_1, 2) + pow(x_2, 2);
  jacobiana[1][1] = 2 * x_1 * x_2 + 2 * x_2 * x_3 - 3 * pow(x_2, 2);
  jacobiana[1][2] = pow(x_2, 2);

  jacobiana[2][0] = 3 * pow(x_1, 2);
  jacobiana[2][1] = 0;
  jacobiana[2][2] = 3 * pow(x_3, 2);

  printf("Valor do chute inicial:[%.2f; %.2f; %.2f]\n", x_1, x_2, x_3);

  printf("Valor para f_1: %f\n", (f_1));
  printf("Valor para f_2: %f\n", (f_2));
  printf("Valor para f_3: %f\n", (f_3));

  printf("======== Matriz Jacobiana ========\n");
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      printf("%.2f ", jacobiana[i][j]);
    }
    printf("\n");
  }

  // para resolver o sistema:

  float mat[3][4];

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      mat[i][j] = jacobiana[i][j];
    }
  }
  mat[0][3] = -(f_1);
  mat[1][3] = -(f_2);
  mat[2][3] = -(f_3);

  printf("======== Matriz Sistema J(xk)S(xk) = -F(xk) ========\n");
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      printf("%.2f ", mat[i][j]);
    }
    printf("\n");
  }

  // Loop do Método de Newton-Raphson
  while (iteracoes < MAX_ITER) {
    EliminacaoGauss(mat, x_1, x_2, x_3);
    
    // Recalcular norma do vetor F para verificar convergência
    f_1 = x_1 * x_3 - x_3 * pow(e, pow(x_1, 2)) + pow(10, -4);
    f_2 = x_1 * (pow(x_1, 2) + pow(x_2, 2)) + pow(x_2, 2) * (x_3 - x_2);
    f_3 = (pow(x_1, 3) + pow(x_3, 3));
    
    norma_f = sqrt(f_1 * f_1 + f_2 * f_2 + f_3 * f_3);
    
    iteracoes++;
    printf("\nIteração %d: Norma F = %.6e\n", iteracoes, norma_f);
    
    if (norma_f < TOL) {
      printf("\nConvergência atingida em %d iterações!\n", iteracoes);
      break;
    }
  }
  
  if (iteracoes >= MAX_ITER) {
    printf("\nNúmero máximo de iterações (%d) atingido sem convergência.\n", MAX_ITER);
  }

  EliminacaoGauss(mat, x_1, x_2, x_3);
  
  // Parar medição de tempo
  fim = clock();
  tempo_processamento = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
  
  printf("\n======== RESUMO FINAL ========\n");
  printf("Total de iterações: %d\n", iteracoes);
  printf("Tempo de processamento: %.6f segundos\n", tempo_processamento);
  printf("Tempo de processamento: %.3f milissegundos\n", tempo_processamento * 1000);
  
  return 0;
}
