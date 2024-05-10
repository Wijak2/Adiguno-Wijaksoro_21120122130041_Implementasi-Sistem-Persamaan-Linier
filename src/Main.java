import java.util.Scanner;
public class Main {
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        System.out.println("Adiguno Wijaksoro");
        System.out.println("21120122130041");
        System.out.println("Program Penyelesaian Persamaan Linier");

        System.out.println("Masukkan jumlah baris dan kolom matriks A:");
        int size = scanner.nextInt();

        double[][] A = new double[size][size];
        double[] b = new double[size];

        System.out.println("Masukkan elemen-elemen matriks A:");
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                A[i][j] = scanner.nextDouble();
            }
        }

        System.out.println("Masukkan elemen vektor b:");
        for (int i = 0; i < size; i++) {
            b[i] = scanner.nextDouble();
        }

        System.out.println("Pilih metode penyelesaian:");
        System.out.println("1. Matriks Balikan");
        System.out.println("2. Lu Gauss");
        System.out.println("3. Crout");
        int choice = scanner.nextInt();

        double[] solution = new double[0];
        switch (choice) {
            case 1:
                solution = penyelesaianMatriksBalikan(A, b);
                break;
            case 2:
                solution = penyelesaianLUGauss(A, b);
                break;
            case 3:
                solution = penyelesaianCrout(A, b);
                break;
            default:
                System.out.println("Pilihan tidak valid.");
                break;
        }

        if (solution == null) {
            System.out.println("Determinan matriks koefisien adalah nol. Persamaan tidak memiliki solusi unik.");
        } else {
            System.out.println("Solusi dari persamaan linear adalah:");
            for (int i = 0; i < size; i++) {
                System.out.println("x" + (i + 1) + " = " + solution[i]);
            }
        }
    }

    // Metode untuk menyelesaikan persamaan linier menggunakan metode Matriks Balikan
    public static double[] penyelesaianMatriksBalikan (double[][] A, double[] b) {
        double[][] A_inverse = inverse(A);
        if (A_inverse == null) {
            return null; // Jika determinan matriks koefisien adalah nol
        }

        // Mengalikan matriks balikan A dengan vektor b
        int size = A.length;
        double[] solution = new double[size];
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                solution[i] += A_inverse[i][j] * b[j];
            }
        }

        return solution;
    }

    // Metode untuk menghitung determinan matriks
    public static double determinant(double[][] A) {
        int n = A.length;
        if (n == 1) {
            return A[0][0];
        } else {
            double det = 0;
            for (int j = 0; j < n; j++) {
                double[][] minor = getMinor(A, 0, j);
                det += Math.pow(-1, 0 + j) * A[0][j] * determinant(minor);
            }
            return det;
        }
    }

    // Metode untuk mendapatkan minor dari suatu matriks
    public static double[][] getMinor(double[][] A, int row, int col) {
        int n = A.length;
        double[][] minor = new double[n - 1][n - 1];
        for (int i = 0, p = 0; i < n; i++) {
            if (i == row) continue;
            for (int j = 0, q = 0; j < n; j++) {
                if (j == col) continue;
                minor[p][q] = A[i][j];
                q++;
            }
            p++;
        }
        return minor;
    }

    // Metode untuk menghitung matriks balikan
    public static double[][] inverse(double[][] A) {
        int n = A.length;
        double[][] A_inverse = new double[n][n];

        // Menghitung determinan matriks A
        double det = determinant(A);
        if (det == 0) {
            return null; // Jika determinan adalah nol, matriks tidak memiliki balikan
        }

        // Menghitung matriks balikan dengan metode adjoin
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double[][] minor = getMinor(A, i, j);
                A_inverse[j][i] = Math.pow(-1, i + j) * determinant(minor) / det;
            }
        }

        return A_inverse;
    }

    // Metode untuk menyelesaikan persamaan linier menggunakan metode LU Gauss
    public static double[] penyelesaianLUGauss(double[][] A, double[] b) {
        int n = A.length;
        double[][] L = new double[n][n];
        double[][] U = new double[n][n];

        // Mendekomposisi matriks A menjadi matriks segitiga atas U dan matriks segitiga bawah L
        if (!luDecomposition(A, L, U)) {
            return null; // Jika matriks koefisien tidak bisa diubah menjadi matriks segitiga atas
        }

        // Menyelesaikan Ly = b menggunakan substitusi maju (forward substitution)
        double[] y = forwardSubstitution(L, b);

        // Menyelesaikan Ux = y menggunakan substitusi mundur (backward substitution)
        return backwardSubstitution(U, y);
    }

    // Metode untuk mendekomposisi matriks A menjadi matriks segitiga atas U dan matriks segitiga bawah L (LU Decomposition)
    public static boolean luDecomposition(double[][] A, double[][] L, double[][] U) {
        int n = A.length;

        for (int i = 0; i < n; i++) {
            for (int k = i; k < n; k++) {
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += (L[i][j] * U[j][k]);
                }
                U[i][k] = A[i][k] - sum;
            }

            if (U[i][i] == 0) {
                return false; // Matriks koefisien tidak bisa diubah menjadi matriks segitiga atas
            }

            for (int k = i; k < n; k++) {
                if (i == k) {
                    L[i][i] = 1;
                } else {
                    double sum = 0;
                    for (int j = 0; j < i; j++) {
                        sum += (L[k][j] * U[j][i]);
                    }
                    L[k][i] = (A[k][i] - sum) / U[i][i];
                }
            }
        }

        return true;
    }

    // Metode untuk menyelesaikan persamaan linier menggunakan metode Crout
    public static double[] penyelesaianCrout(double[][] A, double[] b) {
        int n = A.length;
        double[][] L = new double[n][n];
        double[][] U = new double[n][n];

        // Mendekomposisi matriks A menjadi matriks segitiga atas U dan matriks segitiga bawah L
        if (!croutDecomposition(A, L, U)) {
            return null; // Jika matriks koefisien tidak bisa diubah menjadi matriks segitiga atas
        }

        // Menyelesaikan Ly = b menggunakan substitusi maju (forward substitution)
        double[] y = forwardSubstitution(L, b);

        // Menyelesaikan Ux = y menggunakan substitusi mundur (backward substitution)
        return backwardSubstitution(U, y);
    }

    // Metode untuk mendekomposisi matriks A menjadi matriks segitiga atas U dan matriks segitiga bawah L (Crout Decomposition)
    public static boolean croutDecomposition(double[][] A, double[][] L, double[][] U) {
        int n = A.length;

        for (int j = 0; j < n; j++) {
            U[0][j] = A[0][j];
        }

        for (int i = 0; i < n; i++) {
            L[i][i] = 1;

            for (int j = i; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += L[i][k] * U[k][j];
                }
                U[i][j] = A[i][j] - sum;
            }

            for (int j = i + 1; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++) {
                    sum += L[j][k] * U[k][i];
                }
                if (U[i][i] == 0) {
                    return false; // Matriks koefisien tidak bisa diubah menjadi matriks segitiga atas
                }
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }

        return true;
    }

    // Metode untuk menyelesaikan Ly = b menggunakan substitusi maju (forward substitution)
    public static double[] forwardSubstitution(double[][] L, double[] b) {
        int n = L.length;
        double[] y = new double[n];

        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * y[j];
            }
            y[i] = (b[i] - sum) / L[i][i];
        }

        return y;
    }

    // Metode untuk menyelesaikan Ux = y menggunakan substitusi mundur (backward substitution)
    public static double[] backwardSubstitution(double[][] U, double[] y) {
        int n = U.length;
        double[] x = new double[n];

        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += U[i][j] * x[j];
            }
            x[i] = (y[i] - sum) / U[i][i];
        }

        return x;
    }
}