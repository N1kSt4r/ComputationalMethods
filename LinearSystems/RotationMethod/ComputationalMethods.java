import java.io.IOException;
import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

public class ComputationalMethods {
    class SLEAresult {
        BigDecimal[] solution;
        BigDecimal[] discrepancy;
        BigDecimal det;
        public SLEAresult(BigDecimal[] solution, BigDecimal[] discrepancy, BigDecimal det) {
            this.solution = solution;
            this.det = det;
            this.discrepancy = discrepancy;
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
    public BigDecimal sqrt(BigDecimal x, int Scale) {
        BigDecimal res = new BigDecimal((int)Math.sqrt(x.doubleValue())).setScale(Scale, RoundingMode.HALF_EVEN);
        BigDecimal k = BigDecimal.ONE.setScale(Scale);
        for (int i = 0; i < Scale; ++i) {
            while (res.add(k).pow(2).compareTo(x) <= 0) {
                res = res.add(k);
            }
            k = k.divide(BigDecimal.TEN, RoundingMode.HALF_EVEN);
        }
        return res;
    }
    public BigDecimal sqrt(BigDecimal x) {
        return sqrt(x, 30);
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
    public BigDecimal[] multiply(BigDecimal[][] a, BigDecimal[] b) {
        return matrixToVector(multiply(a, transposition(vectorToMatrix(b))));
    }
    public BigDecimal[] multiply(BigDecimal[] a, BigDecimal[][] b) {
        return matrixToVector(multiply(vectorToMatrix(a), b));
    }
    public BigDecimal[][] vectorToMatrix(BigDecimal[] b) {
        BigDecimal[][] B = new BigDecimal[1][b.length];
        B[0] = Arrays.copyOf(b, b.length);
        return B;
    }
    public BigDecimal[] matrixToVector(BigDecimal[][] b) {
        if (b[0].length != 1) {
            return b[0];
        } else {
            return transposition(b)[0];
        }
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
    public BigDecimal[] solve(BigDecimal[][] matrix, BigDecimal[] vector) {
        BigDecimal[] answer = Arrays.copyOf(vector, vector.length);
        for (int i = matrix.length - 1; i >= 0; --i) {
            answer[i] = answer[i].setScale(30, RoundingMode.HALF_EVEN);
            for (int j = matrix[0].length - 1; j > i; --j) {
                answer[i] = answer[i].subtract(answer[j].multiply(matrix[i][j]));
            }
            answer[i] = answer[i].divide(matrix[i][i], RoundingMode.HALF_EVEN);
        }
        return answer;
    }
    BigDecimal[] eVector(int pos, int n) {
        BigDecimal[] temp = new BigDecimal[n];
        for (int i = 0; i < n; ++i) {
            temp[i] = i == pos ? BigDecimal.ONE : BigDecimal.ZERO;
        }
        return temp;
    }
    BigDecimal[][] eMatrix(int n) {
        BigDecimal[][] temp = new BigDecimal[n][n];
        for (int i = 0; i < n; ++i) {
            temp[i] = eVector(i, n);
        }
        return temp;
    }
    BigDecimal[][] GivensMatrix(BigDecimal cos, BigDecimal sin, int i, int j, int n) {
        BigDecimal[][] temp = eMatrix(n);
        temp[i][i] = cos;
        temp[j][j] = cos;
        temp[i][j] = sin;
        temp[j][i] = sin.negate();
        return temp;
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
    public SLEAresult MethodOfRolling(BigDecimal[][] Matrix, BigDecimal[] Vector) {
        BigDecimal[][] matrix = new BigDecimal[Matrix.length][Matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            matrix[i] = Arrays.copyOf(Matrix[i], Matrix[i].length);
        }
        BigDecimal[] vector = Arrays.copyOf(Vector, Vector.length);

        for (int i = 0; i < matrix.length; ++i) {
            for (int j = i + 1; j < matrix.length; ++j) {
                BigDecimal znam = matrix[i][i].pow(2).add(matrix[j][i].pow(2));
                znam = sqrt(znam);
                BigDecimal cos = matrix[i][i].divide(znam, RoundingMode.HALF_EVEN);
                BigDecimal sin = matrix[j][i].divide(znam, RoundingMode.HALF_EVEN);
                BigDecimal[][] Givens = GivensMatrix(cos, sin, i, j, matrix.length);
                matrix = multiply(Givens, matrix);
                vector = multiply(Givens, vector);
            }
        }
        BigDecimal[] answer = solve(matrix, vector);

        BigDecimal det = BigDecimal.ONE.setScale(30);
        for (int i = 0; i < matrix.length; ++i) {
            det = det.multiply(matrix[i][i]);
        }
        return new SLEAresult(answer, disperancy(Matrix, Vector, answer), det);
    }
}