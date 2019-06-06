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
                System.out.format("%10f  ", elem.setScale(15, RoundingMode.HALF_EVEN).doubleValue());
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
        ComputationalMethods methods = new ComputationalMethods();
        BigDecimal[][] matrixOfSystem = new BigDecimal[4][4];
        BigDecimal[] vector = new BigDecimal[4];

        read(matrixOfSystem, vector, "input.txt");

        methods.methodRelaxing(matrixOfSystem, vector, new BigDecimal("1.07"));
    }
}