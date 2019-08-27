import JAMAJniGsl.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



public final class JAMAJniGslExamples {
	private JAMAJniGslExamples() {}
	public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Exaples for JAMAJniGsl   ### \n \n");
        System.out.println("###   Parameter Preparation   ###");
        
        int M=3, N=3;
        double alpha = 2;
        int matrix_layout = Matrix.LAYOUT.RowMajor;
        long[] pivot;

        int trans = Matrix.TRANSPOSE.Trans;
    	int notrans = Matrix.TRANSPOSE.NoTrans;

        double[][] a = new double[][] {{12,-51,4},{6,167,-68},{-4,24,-41}};
        double[] a1 = new double[] {12, 6,-4};
        double[][] b = new double[][] {{4,-2,-6},{-2,10,9},{-6,9,14}};
        double[] b1 = new double[] {4,-2,-6};
        double[][] c;
        double[][] aa = new double[][] {{1,2,3},{2,4,4}};
        double[][] bb = new double[][] {{1,2},{3,1},{-1,0}};
        //
        //Construct Matrix
        //
        Matrix A = new Matrix(a);
        Matrix B = new Matrix(b);
        Matrix B1 = new Matrix(b1, N);
        Matrix A1 = new Matrix(a1, N);
        Matrix AA = new Matrix(aa);
        Matrix BB = new Matrix(bb);

        Matrix C = new Matrix(M,N);
        Matrix D = new Matrix(M,N);
        Matrix X = new Matrix(M,N);
        Matrix Q = new Matrix(M,N);
        Matrix R = new Matrix(M,N);
        Matrix L = new Matrix(M,N);
        Matrix U = new Matrix(M,N);
        Matrix V = new Matrix(M,N);
        Matrix S = new Matrix(M,N);
		//
        // Print the parameters
        //
        printMatrix("Matrix A", matrix_layout, A.getArray(), M, N);
        printMatrix("Matrix B", matrix_layout, B.getArray(), M, N);
        printMatrix("Matrix B1", matrix_layout, B1.getArray(), N, 1);
        printMatrix("Matrix A1", matrix_layout, A1.getArray(), N, 1);
        System.out.println("M =" + M);
        System.out.println("N =" + N);
        System.out.println("alpha =" + alpha);
        //
        /* ---- Basic matrix algebra ---- */
        //
        System.out.println("\n\n###   Basic matrix algebra   ###");
        System.out.println();
        //
        //Matrix Addition
        //
        System.out.println("\n##  Addition: C = A + B  ##");
        C = A.plus(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //MINUS
        //
        System.out.println("\n##  Subtraction: C = A - B  ##");
        C = A.minus(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //TIMES
        //
        System.out.println("\n##  Multiplication: C = A * B  ##");
        C = A.times(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //TIMES with transpose option
        //
        System.out.println("\n##  Multiplication with transpose option: C = A' * B  ##");
        C = A.times(B,trans,notrans);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //TIMES with transpose option
        //
        System.out.println("\n##  Multiplication with transpose option: C = A * B  ##");
        C = A.times(B,notrans,notrans);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
		//
        //TIMES with transpose option
        //
		printMatrix("Matrix AA", matrix_layout, AA.getArray(), 2, N);
        printMatrix("Matrix BB", matrix_layout, BB.getArray(), M, 2);
        System.out.println("\n##  Multiplication with transpose option: C = AA * BB  ##");
        C = AA.times(BB,notrans,notrans);
        printMatrix("C = ", matrix_layout, C.getArray(), C.getRowDimension(), C.getColumnDimension());
        //
        //TIMES with transpose option
        //
        System.out.println("\n##  Multiplication : C = AA' * BB'  ##");
        C = AA.times(BB,trans,trans);
        printMatrix("C = ", matrix_layout, C.getArray(), C.getRowDimension(), C.getColumnDimension());

	double[][] XXX = {{1,2,3},{0,1,2}};
	Matrix XXXX = new Matrix(XXX);
	Matrix XXXXX = XXXX.times(XXXX,trans,notrans);
	//Matrix XXXX = new Matrix(M,N);
	//XXXX = AA.transpose().times(BB.transpose());
        printMatrix("C = ", matrix_layout, XXXXX.getArray(), XXXXX.getRowDimension(), XXXXX.getColumnDimension());
	//System.out.println(" OK");
        //
        //SCALAR
        //
        System.out.println("\n##  Scalar multiplication: C = alpha * A  ##");
        C = A.times(alpha);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        /* ---- Factorizations ---- */
        //
        System.out.println("\n\n###   Factorizations & Solving Linear Systems   ###");
        System.out.println();
        //
        //LU Decomposition
        //
        System.out.println("\n##  LU Decomposition: A = L*U  ##");
        LUDecomposition LU = A.lu();
        L = LU.getL();
        U = LU.getU();
        pivot = LU.getPivot();
        X = LU.solve(B1);
        printMatrix("L = ", matrix_layout, L.getArray(), M, N);
        printMatrix("U = ", matrix_layout, U.getArray(), M, N);
        printLongArray("The permutation vector pivot is : ", pivot, pivot.length);
        printMatrix("The solution matrix of AX=B1 by LU is", matrix_layout, X.getArray(), N, 1);
        //
        //Cholesky Decomposition
        //
        System.out.println("\n##  Cholesky Decomposition of Matrix B  ##");
        CholeskyDecomposition Chol = B.chol(); 
        L = Chol.getL();
        X = Chol.solve(A1);
        printMatrix("L = ", matrix_layout, L.getArray(), M, N);
        printMatrix("The solution matrix of BX=A1 by Cholesky is", matrix_layout, X.getArray(), N, 1);
        //
        //QR Decomposition
        //
        System.out.println("\n##  QR Decomposition: A = Q*R  ##");
        QRDecomposition QR = A.qr();
        Q = QR.getQ();
        R = QR.getR();
        X = QR.solve(B1);
        printMatrix("Q = ", matrix_layout, Q.getArray(), M, N);
        printMatrix("R = ", matrix_layout, R.getArray(), M, N);
        printMatrix("The solution matrix of AX=B1 by QR is", matrix_layout, X.getArray(), N, 1);
        Matrix test = Q.times(R);
        printMatrix("QR = ", matrix_layout, test.getArray(), M, N);
        //
        //SVD Decomposition
        //
        System.out.println("\n##  SVD Decomposition: A = U*S*V  ##");
        double[] Singtest = {11,22,31,45,52,61,34,24};
        Matrix Sing = new Matrix(Singtest, 4, 2);
        SingularValueDecomposition SVD = Sing.svd();
        U = SVD.getU();
        V = SVD.getV();
        S = SVD.getS();
        Matrix testSVD = Sing.minus(U.times(S.times(V.transpose())));
        printMatrix("U = ", matrix_layout, U.getArray(), 4, 2);
        printMatrix("V = ", matrix_layout, V.getArray(), 2, 2);
        printMatrix("S = ", matrix_layout, S.getArray(), 2, 2);
        printMatrix("test = ", matrix_layout, testSVD.getArray(), 4, 2);

        //
        //Eigenvalue Decomposition
        //
        System.out.println("##  Eigenvalue Decomposition  ##");
        //
        //Asymmetric
        //
        System.out.println("\n#  Asymmetric: A*V = V*D  #");
        EigenvalueDecomposition Eig = A.eig();
        V = Eig.getV();
        D = Eig.getD();
        Matrix Asymme = V.times(D).minus(A.times(V));
        printMatrix("V = ", matrix_layout, V.getArray(), M, N);
        printMatrix("D = ", matrix_layout, D.getArray(), M, N);
        printMatrix("minus = ", matrix_layout, Asymme.getArray(), M, N);
        //
        //Symmetric
        //
        System.out.println("\n#  Symmetric: B = V*D*V'  #");
        Eig = B.eig();
        V = Eig.getV();
        D = Eig.getD();
        Matrix Symme = V.times(D).minus(B.times(V));
        printMatrix("V = ", matrix_layout, V.getArray(), M, N);
        printMatrix("D = ", matrix_layout, D.getArray(), M, N);
        printMatrix("minus = ", matrix_layout, Symme.getArray(), M, N);

        //printMatrix("Matrix A", matrix_layout, A.getArray(), M, N);
        //printMatrix("Matrix B", matrix_layout, B.getArray(), M, N);
        //C = A.plus(B,1);
        //printMatrix("Matrix C", matrix_layout, C.getArray(), M, N);
        //D= C.minus(B.times(1));
        //printMatrix("Matrix D", matrix_layout, D.getArray(), M, N);
    }
    //
    /* Print the matrix X */
    //
    private static void printMatrix(String prompt, int layout, double[][] X, int I, int J) {
        System.out.println(prompt);
        if (layout == Matrix.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j][i]));
                System.out.println();
            }
        }
        else if (layout == Matrix.LAYOUT.RowMajor){
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[i][j]));
                System.out.println();
            }
        }
        else{System.out.println("** Illegal layout setting");}
    }
    //
    /* Print the array X */
    //
    private static void printIntArray(String prompt, int[] X, int L) {
        System.out.println(prompt);
        for (int i=0; i<L; i++) {
            System.out.print("\t" + string(X[i]));
        }
        System.out.println();
    }
    private static void printLongArray(String prompt, long[] X, int L) {
        System.out.println(prompt);
        for (int i=0; i<L; i++) {
            System.out.println("\t" + string(X[i]));
        }
        System.out.println();
    }
    //
    /* Shorter string for real number */
    //
    private static String string(double re) {
        String s="";
        if (re == (long)re)
            s += (long)re;
        else
            s += re;
        return s;
    }
}
