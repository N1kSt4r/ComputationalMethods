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
    public void swap(BigDecimal a[], int i, int j) {
        BigDecimal temp = a[i];
        a[i] = a[j];
        a[j] = temp;
    }
    public void swap(int a[], int i, int j) {
        int temp = a[i];
        a[i] = a[j];
        a[j] = temp;
    }
    public void swapRows(BigDecimal[][] matrix, int i, int j) {
        BigDecimal[] temp = matrix[i];
        matrix[i] = matrix[j];
        matrix[j] = temp;
    }
    public void swapColumns(BigDecimal[][] matrix, int i, int j) {
        BigDecimal temp;
        for (int k = 0; k < matrix.length; ++k) {
           swap(matrix[k], i, j);
        }
    }
    public Point maxOfMatrix(BigDecimal[][] matrix, int from) {
        int row = from;
        int column = from;
        BigDecimal max = matrix[from][from].abs();
        for (int i = from; i < matrix.length; ++i) {
            for (int j = from; j < matrix[0].length; ++j) {
                if (matrix[i][j].abs().compareTo(max) == 1) {
                    max = matrix[i][j].abs();
                    row = i;
                    column = j;
                }
            }
        }
        return new Point(row, column);
    }
    public void sub(BigDecimal[] b, BigDecimal[] a, BigDecimal coef) {
        for (int i = 0; i < b.length; ++i) {
            b[i] = b[i].subtract(a[i].multiply(coef));
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
    public SLEAresult methodGauss(BigDecimal[][] Matrix, BigDecimal[] Vector) {
        BigDecimal[][] matrix = new BigDecimal[Matrix.length][Matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            matrix[i] = Arrays.copyOf(Matrix[i], Matrix[i].length);
        }
        BigDecimal[] vector = Arrays.copyOf(Vector, Vector.length);
        BigDecimal[][] eMatrix = new BigDecimal[matrix.length][matrix.length];
        for (int i = 0; i < eMatrix.length; ++i) {
            for (int j = 0; j < eMatrix.length; ++j) {
                eMatrix[i][j] = i == j ? BigDecimal.ONE : BigDecimal.ZERO;
            }
        }
        int[] columnsSwap = new int[matrix[0].length];
        int[] rowsSwap = new int[matrix.length];
        for (int i = 0; i < columnsSwap.length; ++i) {
            columnsSwap[i] = i;
            rowsSwap[i] = i;
        }

        Point max;
        for (int i = 0; i < matrix.length; ++i) {
            max = maxOfMatrix(matrix, i);
            swapColumns(eMatrix, i, max.column);
            swapRows(eMatrix, i, max.row);
            swapColumns(matrix, i, max.column);
            swapRows(matrix, i, max.row);
            swap(columnsSwap, i, max.column);
            swap(rowsSwap, i, max.row);
            swap(vector, i, max.row);
            for (int j = i + 1; j < matrix.length; ++j) {
                BigDecimal coef = matrix[j][i].divide(matrix[i][i], RoundingMode.HALF_EVEN);
                sub(matrix[j], matrix[i], coef);
                sub(eMatrix[j], eMatrix[i], coef);
                vector[j] = vector[j].subtract(vector[i].multiply(coef));
            }

        }
        int inversions = 0;
        for (int i = 0; i < rowsSwap.length; ++i) {
            for (int j = i + 1; j < rowsSwap.length; ++j) {
                if (rowsSwap[i] > rowsSwap[j]) ++inversions;
                if (columnsSwap[i] > columnsSwap[j]) ++inversions;
            }
        }

        BigDecimal det = BigDecimal.ONE;
        if (inversions % 2 == 1) {
            det = det.negate();
        }

        for (int i = 0; i < matrix.length; ++i) {
            det = det.multiply(matrix[i][i]);
        }

        BigDecimal[] answer = solve(matrix, vector);

        eMatrix = transposition(eMatrix);
        BigDecimal[][] rMatrix = new BigDecimal[eMatrix.length][eMatrix.length];
        for (int i = 0; i < rMatrix.length; ++i) {
            rMatrix[i] = solve(matrix, eMatrix[i]);
        }

        rMatrix = transposition(rMatrix);

        for (int i = 0; i < columnsSwap.length; ++i) {
            swapColumns(rMatrix, i, columnsSwap[i]);
            swap(answer, i, columnsSwap[i]);
            swap(columnsSwap, columnsSwap[i], i);
            for (int j : columnsSwap) {
                System.out.print(j + "  ");
            }
        }
        for (int i = 0; i < rowsSwap.length; ++i) {
            swapRows(rMatrix, i, rowsSwap[i]);
            swap(rowsSwap, rowsSwap[i], i);
        }

        BigDecimal[] discrepancy = Arrays.copyOf(Vector, Vector.length);
        for (int i = 0; i < discrepancy.length; ++i) {
            for (int j = 0; j < Matrix[0].length; ++j) {
                discrepancy[i] = discrepancy[i].subtract(Matrix[i][j].multiply(answer[j]));
            }
        }

        return new SLEAresult(rMatrix, answer, discrepancy, det);
    }
}
