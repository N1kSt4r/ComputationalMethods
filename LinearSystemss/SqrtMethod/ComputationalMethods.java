import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

public class ComputationalMethods {
    class Point {
        public int row;
        public int column;
        Point(int row, int column) {
            this.row = row;
            this.column = column;
        }
    }
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
        BigDecimal[][] c = new BigDecimal[b.length][1];
        for (int i = 0; i < b.length; ++i) {
            c[i][0] = b[i];
        }
        c = multiply(a, c);
        return transposition(c)[0];
    }
    public void swap(BigDecimal a[], int i, int j) {
        BigDecimal temp = a[i];
        a[i] = a[j];
        a[j] = temp;
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
    public BigDecimal[] solveS(BigDecimal[][] S, BigDecimal[] vector) {
        BigDecimal[] y = Arrays.copyOf(vector, vector.length);
        for (int i = 0; i < y.length; ++i) {
            for (int k = 0; k < i; ++k) {
                y[i] = y[i].subtract(S[k][i].multiply(y[k]));
            }
            y[i] = y[i].divide(S[i][i], RoundingMode.HALF_EVEN);
        }

        BigDecimal[] x = Arrays.copyOf(y, y.length);
        for (int i = x.length - 1; i >= 0; --i) {
            for (int k = i + 1; k < x.length; ++k) {
                x[i] = x[i].subtract(S[i][k].multiply(x[k]));
            }
            x[i] = x[i].divide(S[i][i], RoundingMode.HALF_EVEN);
        }
        return x;
    }
    public SLEAresult methodSQRT(BigDecimal[][] Matrix, BigDecimal[] Vector) {
        BigDecimal[][] matrix = multiply(transposition(Matrix), Matrix);
        BigDecimal[] vector = multiply(transposition(Matrix), Vector);
        BigDecimal[][] rMatrix = new BigDecimal[matrix.length][matrix.length];
        for (int i = 0; i < rMatrix.length; ++i) {
            for (int j = 0; j < rMatrix.length; ++j) {
                rMatrix[i][j] = i == j ? BigDecimal.ONE : BigDecimal.ZERO;
            }
            rMatrix[i] = multiply(transposition(Matrix), rMatrix[i]);
        }
        BigDecimal[][] S = new BigDecimal[matrix.length][matrix.length];

        for (int i = 0; i < S.length; ++i) {
            S[i][i] = matrix[i][i];
            for (int k  = 0; k < i; ++k) {
                S[i][i] = S[i][i].subtract(S[k][i].pow(2));
            }
            S[i][i] = sqrt(S[i][i]);

            for (int j = 0; j < i; ++j) {
                S[i][j] = BigDecimal.ZERO;
            }

            for (int j = i + 1; j < S.length; ++j) {
                S[i][j] = matrix[i][j];
                for (int k = 0; k < i; ++k) {
                    S[i][j] = S[i][j].subtract(S[k][j].multiply(S[k][i]));
                }
                S[i][j] = S[i][j].divide(S[i][i], RoundingMode.HALF_EVEN);
            }
        }

        BigDecimal det = BigDecimal.ONE;
        for (int i = 0; i < matrix.length; ++i) {
            det = det.multiply(S[i][i]);
        }

        for (int i = 0; i < rMatrix.length; ++i) {
            rMatrix[i] = solveS(S, rMatrix[i]);
        }
        rMatrix = transposition(rMatrix);

        BigDecimal[] discrepancy = Arrays.copyOf(vector, vector.length);
        vector = solveS(S, vector);
        for (int i = 0; i < discrepancy.length; ++i) {
            for (int j = 0; j < Matrix[0].length; ++j) {
                discrepancy[i] = discrepancy[i].subtract(matrix[i][j].multiply(vector[j]));
            }
        }

        return new SLEAresult(rMatrix, vector, discrepancy, det);
    }
}