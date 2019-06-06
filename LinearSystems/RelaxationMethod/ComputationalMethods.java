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
    public void iterationGayssZeidel(BigDecimal[] answer, BigDecimal[][] B, BigDecimal[] g, BigDecimal w) {
        for (int i = 0; i < answer.length; ++i) {
            BigDecimal sum = BigDecimal.ZERO;
            for (int j = 0; j < answer.length; ++j) {
                sum = sum.add(B[i][j].multiply(answer[j]));
            }
            answer[i] = answer[i].multiply(BigDecimal.ONE.subtract(w));
            answer[i] = answer[i].add(sum.multiply(w));
            answer[i] = answer[i].add(g[i].multiply(w));
        }
    }
    public void methodRelaxing(BigDecimal[][] Matrix, BigDecimal[] Vector, BigDecimal w) {
        BigDecimal[][] matrix = new BigDecimal[Matrix.length][Matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            matrix[i] = Arrays.copyOf(Matrix[i], Matrix[i].length);
        }
        BigDecimal[] answerRelax = Arrays.copyOf(Vector, Vector.length);

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
        final BigDecimal EPS = new BigDecimal("0.1").pow(29);

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

        int countIterationRelax = 0;
        while (rate(disperancy(Matrix, Vector, answerRelax)).compareTo(EPS) > 0) {
            iterationGayssZeidel(answerRelax, B, g, w);
            ++countIterationRelax;
        }
        System.out.println("Реальное количество итераций в методе Релаксации: " + countIterationRelax + '\n');

        System.out.print("Вектор решений: ");
        Main.display(answerRelax);
        System.out.print("Вектор невязки: ");
        Main.displayE(disperancy(Matrix, Vector, answerRelax));
    }
}