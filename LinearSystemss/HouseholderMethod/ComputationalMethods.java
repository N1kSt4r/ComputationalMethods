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
    BigDecimal[] vectorSub(BigDecimal[] a, BigDecimal[] b, BigDecimal coef) {
        BigDecimal[] temp = new BigDecimal[a.length];
        for (int i = 0; i < a.length; ++i) {
            temp[i] = a[i].subtract(b[i].multiply(coef));
        }
        return temp;
    }
    BigDecimal[][] matrixSub(BigDecimal[][] a, BigDecimal[][] b, BigDecimal coef) {
        BigDecimal[][] temp = new BigDecimal[a.length][];
        for (int i = 0; i < a.length; ++i) {
            temp[i] = vectorSub(a[i], b[i], coef);
        }
        return temp;
    }
    BigDecimal[] multiply(BigDecimal[] a, BigDecimal b) {
        BigDecimal[] temp = new BigDecimal[a.length];
        for (int i = 0; i < a.length; ++i) {
            temp[i] = a[i].multiply(b);
        }
        return temp;
    }
    BigDecimal scalarMultiply(BigDecimal[] a, BigDecimal[] b) {
        BigDecimal sum = BigDecimal.ZERO;
        for (int i = 0; i < a.length; ++i) {
            sum = sum.add(a[i].multiply(b[i]));
        }
        return sum;
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
    public SLEAresult Householder(BigDecimal[][] Matrix, BigDecimal[] Vector) {
        BigDecimal[][] matrix = new BigDecimal[Matrix.length][Matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            matrix[i] = Arrays.copyOf(Matrix[i], Matrix[i].length);
        }
        BigDecimal[] vector = Arrays.copyOf(Vector, Vector.length);
        for (int i = 0; i < matrix.length - 1; ++i) {
            BigDecimal[] s = transposition(matrix)[i];
            for (int j = 0; j < i; ++j) {
                s[j] = BigDecimal.ZERO;
            }
            BigDecimal afla = sqrt(scalarMultiply(s, s), 30);
            BigDecimal[] sMinusAlfaE = vectorSub(s, multiply(eVector(i, vector.length), afla), BigDecimal.ONE);

            BigDecimal x = BigDecimal.ONE.setScale(30).divide(
                    sqrt(scalarMultiply(s, sMinusAlfaE).multiply(new BigDecimal("2")), 30), RoundingMode.HALF_EVEN);

            BigDecimal[] w = multiply(sMinusAlfaE, x);
            BigDecimal[][] V = matrixSub(eMatrix(matrix.length),
                    multiply(transposition(vectorToMatrix(w)), vectorToMatrix(w)), new BigDecimal("2"));

            matrix = multiply(V, matrix);

            vector = multiply(V, vector);
        }
        BigDecimal det = BigDecimal.ONE.setScale(30);
        for (int i = 0; i < matrix.length; ++i) {
            det = det.multiply(matrix[i][i]);
        }
        return new SLEAresult(solve(matrix, vector), disperancy(Matrix, Vector, solve(matrix, vector)), det.abs());
    }
}