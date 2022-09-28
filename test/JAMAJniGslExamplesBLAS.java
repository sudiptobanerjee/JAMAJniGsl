import JAMAJniGsl.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;


public final class JAMAJniGslExamplesBLAS {
	private JAMAJniGslExamplesBLAS() {}
	public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
		int M=2, N=4, K=3;
		float alpha = 2, beta = -1;
		double alphad = 2, betad = -1;

		int LAYOUT = Matrix.LAYOUT.ColMajor;    
		int transA = Matrix.TRANSPOSE.NoTrans;
		int transB = Matrix.TRANSPOSE.NoTrans;
		int uplo;
		int diag;
		int side;

		double resultd;

		double[] dA = new double[] {1,2,3,4,5,6};
		double[] cpdA = new double[] {1,5,4,3,2,6};
        double[] dotd = new double[] {1,0,1,-1};//with dim M
        double[] cpdotd = new double[] {5,5,13,9};//with dim M
        double[] dotr = new double[] {1,0};
        double[] cpdotr = new double[] {1,0};
        double[] doti = new double[] {1,-1};
        double[] cpdoti = new double[] {1,-1};
        double[] xd = new double[] {1,2,3};
        double[] cpxd = new double[] {1,2,3};
        double[] yd = new double[] {0,1,-1};
        double[] cpyd = new double[] {0,1,-1};
        double[] trid = new double[] {1,2,3,0,4,5,0,0,0,};
        double[] solvetrid = new double[] {1,2,3,0,4,5,0,0,1};
        double[] symfulld = new double[]{1,4,0,4,2,5,0,5,3};
        double[] Bd = new double[] {7,8,9,10,11,12,13,14,15,16,17,18};//with dim K*N
        double[] Cd = new double[] {2,-4,3,7,-3,-7,0,18};//with dim M*N
        double[] myad = new double[] {1,4,2,5,3,6};//dim M*K
        double[] cpmyad = new double[] {1,4,2,5,3,6};//dim M*K
   


        /* Level 1: */
        System.out.println("Level 1:");
	    //dscal
        System.out.println();
        System.out.println("dscal: scales a vector by a constant");
        CprintMatrix("vector x = ", dA, M*K, 1);
        System.out.println("alpha=" + string(alpha));

	    //N = M*K
	    //alpha = alphad	
	    //x = dA

        Matrix.dscal(M*K, alphad, dA);
        CprintMatrix("Afterscaling: vector x = ", dA, M*K, 1);

        Matrix.dscal(M*K, 1/alphad, dA);

        //daxpy
        System.out.println();
        System.out.println("daxpy: y = a * x + y ");
        System.out.println("a = " + string(alphad));
        CprintMatrix("x = ", cpdotr, M, 1);
        doti = cpdoti.clone();
        CprintMatrix("y = ", doti, M, 1);

        //N = M
        //alpha = alphad
        //X = cpdotr
        //incX = 1
        //Y = doti
        //incY = 1

        Matrix.daxpy(M, alphad, cpdotr, doti);
        CprintMatrix("Result = ", doti, M, 1);

        // ddot
        System.out.println();
        System.out.println("ddot: ddot = inner product of vector x and y");
        CprintMatrix("vector x = ", dA, M*K, 1);
        CprintMatrix("vector y = ", cpdA, M*K, 1);

        //N = M*k
        //X = dA
        //incX = 1
        //Y = cpdA
        //incY = 1

        resultd = Matrix.ddot(M*K, dA, cpdA);
        System.out.printf("result = %f\n", resultd);

        /* Level 2: */
        System.out.println();
        System.out.printf("Level 2:");
        System.out.println();
        //dgemv
        System.out.println();
        System.out.println("dgemv: y = alpha * A * x + beta * y , or y = alpha * (A**T) * x + beta * y");
        System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
        dotd = cpdotd.clone();
        doti = cpdoti.clone();
        dotr = cpdotr.clone();
        RprintMatrix("A = ", dotd, M, M);
        RprintMatrix("x = ", dotr, M, 1);
        RprintMatrix("y = ", doti, M, 1);
        LAYOUT = Matrix.LAYOUT.RowMajor;    

        //LAYOUT = LAYOUT
        //TransA = transA
        //M = M
        //N = M
        //alpha = alpha
        //A = dotd
        //lda = M
        //x = dotr
        //beta = beta
        //y = doti
        //incY = 1

        transA = Matrix.TRANSPOSE.Trans;
        Matrix.dgemv(transA, M, M, alpha, dotd, dotr, beta, doti);
        RprintMatrix("Result: y = ", doti, M, 1);

        //dtrmv
        System.out.println();
        System.out.println("dtrmv: x = A * x or x = A**T * x");
        RprintMatrix("A = ", trid, K, K);
        RprintMatrix("x = ", xd, K, 1);
        LAYOUT = Matrix.LAYOUT.RowMajor;    
        transA = Matrix.TRANSPOSE.NoTrans;
        uplo = Matrix.UPLO.Upper;
        diag = Matrix.DIAG.NonUnit;

        //LAYOUT = LAYOUT
        //Uplo = uplo
        //TransA = transA
        //Diag = diag
        //N = K 
        //A = trid
        //lda = K
        //x = xd
        //incX = 1

        Matrix.dtrmv( uplo, transA, diag, K, trid, xd);
        RprintMatrix("Result: x = ", xd, K, 1);

        //dsymv
        System.out.println();
        System.out.println("dsymv: y = alpha * A * x + beta * y, where A is symmetric matrix");
        System.out.printf("alpha = %f\n", alpha);
        System.out.printf("beta = %f\n", beta);
        RprintMatrix("A = ", symfulld, K, K);
        xd = cpxd.clone();
        yd = cpyd.clone();
        RprintMatrix("x = ", xd, K, 1);
        RprintMatrix("y = ", yd, K, 1);
        LAYOUT = Matrix.LAYOUT.RowMajor;    
        uplo = Matrix.UPLO.Upper;

        //uplo = uplo
        //N = K
        //alpha = alpha
        //a = symfulld
        //lda = K
        //x = xd
        //incx = 1
        //beta = beta
        //y = yd
        //incy = 1

        Matrix.dsymv(uplo, K, alpha, symfulld, xd, beta, yd);
        RprintMatrix("Result: x = ", yd, K, 1);
        

        /* Level 3: */
        System.out.println();
        System.out.printf("Level 3:");
        System.out.println();

        //dgemm
        System.out.println();
        System.out.println("dgemm: C = alpha*op( A )*op( B ) + beta*C ");
        System.out.println("transA = N, transB = N, LAYOUT = C");
        System.out.println("alpha = " + string(alpha));
        System.out.println("beta = " + string(beta));
        Matrix.dscal(M*K, 2, dA);
        RprintMatrix("Matrix A", dA, M, K);
        RprintMatrix("Matrix B", Bd, K, N);
        RprintMatrix("Matrix C", Cd, M, N);
        transA = Matrix.TRANSPOSE.NoTrans;
        transB = Matrix.TRANSPOSE.NoTrans;
        LAYOUT = Matrix.LAYOUT.RowMajor;

        //TransA = transA
        //TransB = transB
        //M = M
        //N = N
        //K = K
        //alpha = alphad    
        //a = dA
        //lda = M
        //b = Bd
        //ldb = K
        //beta = betad
        //c = Cd
        //ldc = M

        Matrix.dgemm(transA, transB, M, K, K, N, alphad, dA, Bd, betad, Cd);
        RprintMatrix("Resulting C", Cd, M, N);

	// test case involving one of the product matrices being transposed
	Cd = new double[6];
        transA = Matrix.TRANSPOSE.Trans;
        Matrix.dgemm(transA, transB, 2, 3, 2, 2, 1.0, dA, dotd, 0.0, Cd);
        RprintMatrix("Resulting C", Cd, 3, 2);


        //dtrmm
        System.out.println();
        System.out.println("dtrmm: B = alpha * op(A) * B, or B = alpha * B * op(A), where A is a triangular matrix");
        System.out.println("alpha = " + string(alpha));
        RprintMatrix("Matrix A", trid, K, K);
        RprintMatrix("Matrix B", myad, K, M);
        LAYOUT = Matrix.LAYOUT.RowMajor;
        transA = Matrix.TRANSPOSE.NoTrans;
        uplo = Matrix.UPLO.Upper;
        side = Matrix.SIDE.Left;
        diag = Matrix.DIAG.NonUnit;

        //Side = side
        //Uplo = uplo
        //Trans = transA
        //Diag = diag
        //M = K
        //N = M
        //alpha = alphad
        //a = trid
        //lda = K
        //b = myad
        //ldb = K

        Matrix.dtrmm(side, uplo, transA, diag, K, M, alphad, trid, myad);
        RprintMatrix("Resulting C", myad, K, M);

        //dtrsm
        System.out.println();
        System.out.println("dtrsm: op(A) * X = alpha * B, or X * op(A) = alpha * B, where A is a triangular matrix");
        System.out.println("alpha = " + string(alpha));
        RprintMatrix("Matrix A", solvetrid, K, K);
        myad = cpmyad.clone();
        RprintMatrix("Matrix B", myad, K, M);
        LAYOUT = Matrix.LAYOUT.RowMajor;
        transA = Matrix.TRANSPOSE.NoTrans;
        uplo = Matrix.UPLO.Upper;
        side = Matrix.SIDE.Left;
        diag = Matrix.DIAG.NonUnit;

        //Side = side
        //Uplo = uplo
        //TransA = transA
        //Diag = diag
        //M = K
        //N = M
        //alpha = alphad
        //a = solvetrid
        //lda = K
        //b = myad
        //ldb = K

        Matrix.dtrsm(side, uplo, transA, diag, K, M, alphad, solvetrid, myad);
        RprintMatrix("Solution X = ", myad, K, M);

        //dsymm
        System.out.println();
        System.out.println("dsymm: C = alpha * A * B + beta * C, or C = alpha * B * A + beta * C ");
        System.out.println("alpha = " + string(alpha));
        System.out.println("beta = " + string(beta));
        RprintMatrix("Matrix A", symfulld, K, K);
        RprintMatrix("Matrix B", cpdA, K, M);
        RprintMatrix("Matrix C", dA, K, M);
        LAYOUT = Matrix.LAYOUT.RowMajor;
        side = Matrix.SIDE.Left;
        uplo = Matrix.UPLO.Upper;

        //side = side
        //uplo = uplo
        //M = K
        //N = M
        //alpha = alphad
        //a = symfulld
        //lda = K
        //b = cpdA
        //ldb = K
        //beta = betad
        //c = dA
        //ldc = K

        Matrix.dsymm(side, uplo, K, M, alphad, symfulld, cpdA, betad, dA);
        RprintMatrix("Resulting C", dA, K, M);
    }

    //print vector into matrix, assume vector is col-majored
    private static void CprintMatrix(String prompt, double[] X, int I, int J) {
    	System.out.println(prompt);
    	double[][] mat = new double[I][J]; 
    	int count = 0;

    	for (int j=0; j<J; j++){
    		for (int i=0; i<I; i++){
    			mat[i][j] = X[count];
    			count = count+1;
    		}
    	}

    	for (int i=0; i<I; i++) {
    		for (int j=0; j<J; j++)
                //System.out.print("\t" + string(mat[i][j]));
    			System.out.printf("%.4f\t", mat[i][j]);
    		System.out.println();
    	}
    }

    private static void CprintMatrix(String prompt, float[] X, int I, int J) {
    	System.out.println(prompt);
    	float[][] mat = new float[I][J]; 
    	int count = 0;

    	for (int j=0; j<J; j++){
    		for (int i=0; i<I; i++){
    			mat[i][j] = X[count];
    			count = count+1;
    		}
    	}

    	for (int i=0; i<I; i++) {
    		for (int j=0; j<J; j++)
                    //System.out.print("\t" + string(mat[i][j]));
    			System.out.printf("%.4f\t", mat[i][j]);
    		System.out.println();
    	}
    }

        //print vector into matrix, assume vector is row-majored
    private static void RprintMatrix(String prompt, double[] X, int I, int J) {
    	System.out.println(prompt);
    	for (int i=0; i<I; i++) {
    		for (int j=0; j<J; j++)
    			System.out.printf("%.4f\t", X[i*J+j]);
    		System.out.println();
    	}
    }

    private static void RprintMatrix(String prompt, float[] X, int I, int J) {
    	System.out.println(prompt);
    	for (int i=0; i<I; i++) {
    		for (int j=0; j<J; j++)
    			System.out.printf("%.4f\t", X[i*J+j]);
    		System.out.println();
    	}
    }

    /** Shorter string for real number. */
    private static String string(double re) {
    	String s="";
    	if (re == (long)re)
    		s += (long)re;
    	else
    		s += re;
    	return s;
    }
}
