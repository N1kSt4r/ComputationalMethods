import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

import static java.lang.Math.floor;
import static java.lang.Math.log;

public class ComputationalMethods {
    BigDecimal rate(BigDecimal[] vector) {
        BigDecimal rate = BigDecimal.ZERO;
        for (BigDecimal elem : vector) {
            rate = rate.add(elem.abs());
        }
        return rate;
    }
    BigDecimal rate(BigDecimal[][] matrix) {
        BigDecimal rate = BigDecimal.ZERO;
        for (int j = 0; j < matrix.length; ++j) {
            BigDecimal temp = rate(matrix[j]);
            rate = rate.compareTo(temp) < 0 ? temp : rate;
        }
        return rate;
    }
    BigDecimal[] disperancy(BigDecimal[][] Matrix, BigDecimal[] Vector, BigDecimal[] answer) {
        BigDecimal[] discrepancy = Arrays.copyOf(Vector, Vector.length);
        for (int i = 0; i < discrepancy.length; ++i) {
            for (int j = 0; j < Matrix[0].length; ++j) {
                discrepancy[i] = discrepancy[i].subtract(Matrix[i][j].multiply(answer[j]));
            }
        }
        return discrepancy;
    }
    public void iterationGayssZeidel(BigDecimal[] answer, BigDecimal[][] B, BigDecimal[] g) {
        for (int i = 0; i < answer.length; ++i) {
            answer[i] = BigDecimal.ZERO;
            for (int j = 0; j < answer.length; ++j) {
                answer[i] = answer[i].add(B[i][j].multiply(answer[j]));
            }
            answer[i] = answer[i].add(g[i]);
        }
    }
    public BigDecimal[] iterationJacobi(BigDecimal[] vector, BigDecimal[][] B, BigDecimal[] g) {
        BigDecimal[] answer = new BigDecimal[vector.length];
        for (int i = 0; i < answer.length; ++i) {
            answer[i] = BigDecimal.ZERO;
            for (int j = 0; j < answer.length; ++j) {
                answer[i] = answer[i].add(B[i][j].multiply(vector[j]));
            }
            answer[i] = answer[i].add(g[i]);
        }
        return answer;
    }
    public void methodJacobi(BigDecimal[][] Matrix, BigDecimal[] Vector) {
        BigDecimal[][] matrix = new BigDecimal[Matrix.length][Matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            matrix[i] = Arrays.copyOf(Matrix[i], Matrix[i].length);
        }
        BigDecimal[] answerJacobi = Arrays.copyOf(Vector, Vector.length);
        BigDecimal[] answerGayssZeidel = Arrays.copyOf(Vector, Vector.length);

        BigDecimal[][] B = new BigDecimal[Matrix.length][Matrix.length];
        for (int i = 0; i < B.length; ++i) {
            for (int j = 0; j < B.length; ++j) {
                if (i != j) {
                    B[i][j] = matrix[i][j].divide(matrix[i][i], RoundingMode.HALF_EVEN).negate();
                } else {
                    B[i][j] = BigDecimal.ZERO;
                }
            }
        }
        BigDecimal[] g = new BigDecimal[Vector.length];
        for (int i = 0; i < g.length; ++i) {
            g[i] = Vector[i].divide(matrix[i][i], RoundingMode.HALF_EVEN);
        }
        final BigDecimal EPS = new BigDecimal("0.1").pow(4);

        System.out.println("Исходная матрица: ");
        Main.display(Matrix);
        System.out.println();
        System.out.println("Столбец свободных членов: ");
        Main.display(Vector);
        System.out.println();

        System.out.println("Матрица В: ");
        Main.display(B);
        System.out.println();
        System.out.println("Вектор g: ");
        Main.display(g);
        System.out.println();

        System.out.println("Норма B: " + rate(B).doubleValue());
        System.out.println("Норма g: " + rate(g).doubleValue() + "\n");

        System.out.format("Предположительное количество итераций: %.0f\n",
                floor((log(EPS.doubleValue() * (1 - rate(B).doubleValue()) / rate(g).doubleValue()))
                        / log(rate(B).doubleValue())));

        int countIterationJacobi = 0;
        while (rate(disperancy(Matrix, Vector, answerJacobi)).compareTo(EPS) > 0) {
            answerJacobi = iterationJacobi(answerJacobi, B, g);
            ++countIterationJacobi;
        }
        System.out.println("Реальное количество итераций в методе Якоби: " + countIterationJacobi);

        int countIterationGayssZeidel = 0;
        while (rate(disperancy(Matrix, Vector, answerGayssZeidel)).compareTo(EPS) > 0) {
            iterationGayssZeidel(answerGayssZeidel, B, g);
            ++countIterationGayssZeidel;
        }
        System.out.println("Реальное количество итераций в методе Гаусса-Зейделя: " + countIterationGayssZeidel + '\n');

        System.out.print("Вектор решений Якоби: ");
        Main.display(answerJacobi);
        System.out.print("Вектор невязки Якоби: ");
        Main.displayE(disperancy(Matrix, Vector, answerJacobi));
        System.out.println();

        System.out.print("Вектор решений Гаусса-Зейделя: ");
        Main.display(answerGayssZeidel);
        System.out.print("Вектор невязки Гаусса-Зейделя: ");
        Main.displayE(disperancy(Matrix, Vector, answerGayssZeidel));
    }
}