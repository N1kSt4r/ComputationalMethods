import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

public class ComputationalMethods {
    class SLEAresult {
        BigDecimal[][] rMatrix;
        BigDecimal[] solution;
        BigDecimal[] discrepancy;
        BigDecimal det;
        public SLEAresult(BigDecimal[][] rMatrix, BigDecimal[] solution, BigDecimal[] discrepancy, BigDecimal det) {
            this.rMatrix = rMatrix;
            this.solution = solution;
            this.det = det;
            this.discrepancy = discrepancy;
        }
        public BigDecimal[][] getInversedMatrix() {
            return rMatrix;
        }
        public BigDecimal[] getSolution(){
            return solution;
        }
        public BigDecimal[] getDiscrepancy() {
            return discrepancy;
        }
        public BigDecimal getDet() {
            return det;
        }
    }
    public BigDecimal[][] multiply(BigDecimal[][] a, BigDecimal[][] b) {
        BigDecimal[][] c = new BigDecimal[a.length][b[0].length];
        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] = BigDecimal.ZERO.setScale(30, RoundingMode.HALF_EVEN);
                for (int k = 0; k < b.length; ++k) {
                    c[i][j] = c[i][j].add(a[i][k].multiply(b[k][j]));
                }
            }
        }
        return c;
    }
    public BigDecimal[][] transposition(BigDecimal[][] a) {
        BigDecimal[][] b = new BigDecimal[a[0].length][a.length];
        for (int i = 0; i < b.length; ++i) {
            for (int j = 0; j < b[0].length; ++j) {
                b[i][j] = a[j][i];
            }
        }
        return b;
    }
    public BigDecimal[] solveSweep(BigDecimal[][] matrix, BigDecimal[] vector) {
        BigDecimal[] answer = Arrays.copyOf(vector, vector.length);
        int n = answer.length - 1;
        BigDecimal[] alfa = new BigDecimal[answer.length - 1];
        BigDecimal[] beta = new BigDecimal[answer.length];
        alfa[0] = matrix[0][1].negate().divide(matrix[0][0], RoundingMode.HALF_EVEN);
        beta[0] = answer[0].divide(matrix[0][0], RoundingMode.HALF_EVEN);
        for (int i = 1; i < n; ++i) {
            BigDecimal znam = matrix[i][i].subtract(
                    alfa[i - 1].multiply(matrix[i][i - 1]).negate());
            alfa[i] = matrix[i][i + 1].negate().divide(znam, RoundingMode.HALF_EVEN);
            beta[i] = (answer[i].add(matrix[i][i - 1].negate().multiply(beta[i - 1])))
                    .divide(znam, RoundingMode.HALF_EVEN);
        }
        beta[n] = (answer[n].add(matrix[n][n - 1].negate().multiply(beta[n - 1])))
                .divide(matrix[n][n].subtract(alfa[n - 1].multiply(matrix[n][n - 1]).negate()),
                        RoundingMode.HALF_EVEN);

        answer[n] = beta[n];
        for (int i = n - 1; i >= 0; --i) {
            answer[i] = beta[i].add(answer[i + 1].multiply(alfa[i]));
        }
        return answer;
    }
    public SLEAresult methodOfSweep(BigDecimal[][] Matrix, BigDecimal[] Vector) {
        BigDecimal[][] matrix = new BigDecimal[Matrix.length][Matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            matrix[i] = Arrays.copyOf(Matrix[i], Matrix[i].length);
        }
        BigDecimal[] answer = Arrays.copyOf(Vector, Vector.length);
        BigDecimal[][] rMatrix = new BigDecimal[matrix.length][matrix.length];
        for (int i = 0; i < rMatrix.length; ++i) {
            for (int j = 0; j < rMatrix.length; ++j) {
                rMatrix[i][j] = i == j ? BigDecimal.ONE : BigDecimal.ZERO;
                rMatrix[i][j] = rMatrix[i][j].setScale(30, RoundingMode.HALF_EVEN);
            }
        }

        int n = answer.length - 1;
        BigDecimal det = matrix[0][0];
        BigDecimal[] alfa = new BigDecimal[answer.length];
        alfa[0] = matrix[0][1].negate().divide(matrix[0][0], RoundingMode.HALF_EVEN);
        for (int i = 1; i < n; ++i) {
            BigDecimal znam = matrix[i][i].subtract(
                    alfa[i - 1].multiply(matrix[i][i - 1]).negate());
            alfa[i] = matrix[i][i + 1].negate().divide(znam, RoundingMode.HALF_EVEN);
            det = det.multiply(znam);
        }
        det = det.multiply(matrix[n][n].subtract(
                alfa[n - 1].multiply(matrix[n][n - 1]).negate()));

        answer = solveSweep(matrix, answer);

        for (int i = 0; i < rMatrix.length; ++i) {
            rMatrix[i] = solveSweep(matrix, rMatrix[i]);
        }
        rMatrix = transposition(rMatrix);

        BigDecimal[] discrepancy = Arrays.copyOf(Vector, Vector.length);
        for (int i = 0; i < discrepancy.length; ++i) {
            for (int j = 0; j < Matrix[0].length; ++j) {
                discrepancy[i] = discrepancy[i].subtract(Matrix[i][j].multiply(answer[j]));
            }
        }

        return new SLEAresult(rMatrix, answer, discrepancy, det);
    }
}