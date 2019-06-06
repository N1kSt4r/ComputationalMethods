import java.io.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Scanner;

public class Main {
    static void read(BigDecimal[][] matrix, BigDecimal[] vector, String fileName) {
        try (Scanner in = new Scanner(new File(fileName))) {
            for (BigDecimal[] row : matrix) {
                for (int i = 0; i < row.length; ++i) {
                    row[i] = new BigDecimal(in.next());
                    row[i] = row[i].setScale(30, RoundingMode.HALF_EVEN);
                }
            }
            for (int i = 0; i < vector.length; ++i) {
                vector[i] = new BigDecimal(in.next());
                vector[i] = vector[i].setScale(30, RoundingMode.HALF_EVEN);
            }
        } catch (IOException err) {
            System.out.println(err.getMessage());
        }
    }
    static public void display(BigDecimal[][] matrix) {
        for (BigDecimal[] row : matrix) {
            for (BigDecimal elem : row) {
                System.out.format("%10.5f  ", elem.setScale(15, RoundingMode.HALF_EVEN).doubleValue());
            }
            System.out.println();
        }
    }
    static public void display(BigDecimal[] vector) {
        for (BigDecimal elem : vector) {
            System.out.format("%-15f", elem.doubleValue());
        }
        System.out.println();
    }
    static public void displayE(BigDecimal[][] matrix) {
        for (BigDecimal[] row : matrix) {
            for (BigDecimal elem : row) {
                System.out.format(" %+10e ", elem.doubleValue());
            }
            System.out.println();
        }
    }
    static public void displayE(BigDecimal[] vector) {
        for (BigDecimal elem : vector) {
            System.out.format("%-15e", elem.doubleValue());
        }
        System.out.println();
    }
    static public void main(String[] args) {
        int n = 4;
        ComputationalMethods methods = new ComputationalMethods();
        BigDecimal[][] matrixOfSystem = new BigDecimal[n][n];
        BigDecimal[] vector = new BigDecimal[n];

        read(matrixOfSystem, vector, "input.txt");

        ComputationalMethods.SLEAresult result =
                methods.Householder(matrixOfSystem, vector);

        System.out.println("Матрица системы:");
        display(matrixOfSystem);
        System.out.println();

        System.out.println("Столбец свободных членов:");
        display(vector);
        System.out.println();

        System.out.println("Вектор решений:");
        display(result.getSolution());
        System.out.println();

        System.out.println("Модуль определителя:\n" + result.getDet().doubleValue() + '\n');

        System.out.println("Невязка:");
        displayE(result.getDiscrepancy());
        System.out.println();
    }
}