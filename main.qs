namespace Least.Squares.Solver {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Simulation;
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Synthesis;
    open Microsoft.Quantum.Preparation;



    //Prep (DONE)
        //Prep b
        //Prep A (Hermitian)
    
    
    //Get U given A an t
        //U = e^iAt (DONE!)
        //Find U_daggar

    //Convert U to Weird Variables
        //Find theta, phi, lambda, and gamma such that
        //Equation 38 in walkthrough
        //For U, U^2, U^-1, and (U^-1)^2 etc for all powers
    
    //Controlled Rotations
        //Find theta for all eigenvalues (2 arcsin(1 / c)) look in walkthrough for equation for c
        //Make Ry(theta, ) gate
        //Controlled Ry(controlbits, (theta, ancilla))
    

    //Main
        //Do all of the above functions in order
        //Uncompute QPE
    
    

    operation prepareStateB (data : Double[][], register : Qubit[]) : Unit {
        //This prepares the vector |b> as the amplitudes of the |1> state of an ancillary qubit
        //when entangled with the index qubit. ie: b[5] = |5>0.7|1>
        //This is the method used in the 25 page paper.

        mutable sumOfSquares = 0.;
        for i in data {set sumOfSquares += i[1] ^ 2.;} //Find largest element of b
        let normFactor = Sqrt(sumOfSquares);

        mutable b = [];
        for i in data {set b += [i[1] / normFactor];} //Gets the b vecotor normalized w/ all entries < 1
        Message($"{b}");

        PrepareArbitraryStateD(b, LittleEndian(register));
    }


    function displayMatrix(inputMatrix : Double[][], name : String) : Unit {
        let matrix = transpose(inputMatrix);
        mutable o = name + "\n|";

        mutable maxLength = 0;
        for column in matrix {
            for element in column {
                if element > 1. and (Ceiling(Log10(element))) > maxLength {set maxLength = Ceiling(Log10(element));}
            }
        }

        for column in matrix {
            for element in column {
                set o += DoubleAsString(element);
                for space in 0 .. maxLength - (element > 1. ? Ceiling(Log10(element)) | 1) {set o += " ";}
            }
            set o += "|\n|";
        }
        Message(o);
    }


    function transpose(matrix : Double[][]) : Double[][] {
        mutable o = [[0.], size=0];
        for i in 0 .. Length(matrix[0]) - 1 {
            mutable temp = [0., size=0];
            for j in 0 .. Length(matrix) - 1 {
                set temp += [matrix[j][i]];
            }
            set o += [temp];
        }
        return o;
    }
    function transposeC(matrix : Complex[][]) : Complex[][] {
        mutable o = [[Complex(0., 0.)], size=0];
        for i in 0 .. Length(matrix[0]) - 1 {
            mutable temp = [Complex(0., 0.), size=0];
            for j in 0 .. Length(matrix) - 1 {
                set temp += [matrix[j][i]];
            }
            set o += [temp];
        }
        return o;
    }


    function prepareOriginalMatrixA (data : Double[][], width : Int) : Double[][] {
        mutable A = [[0.], size=0];
        for column in 0 .. width - 1 {
            mutable temp = [0., size=0];
            for row in data {
                set temp += [row[0] ^ IntAsDouble(width - column - 1)];
            }
            set A += [temp];
        }
        return A;
    }


    function convertAtoHermitian (A : Double[][]) : Double[][] {
        //Converts A to
        //(0  A)
        //(At 0)

        let At = transpose(A); //Should technically be conjugate transpose, but it is real-valued
        let (w, h) = (Length(A), Length(A[0]));
        mutable o = [[0.], size=0];

        for i in 0 .. w + h - 1 {
            mutable temp = [0., size=0];
            for j in 0 .. h + w - 1 {
                if i >= h and j < h {
                    set temp += [A[i - h][j]];
                }
                elif i < h and j >= h {
                    set temp += [At[i][j - h]];
                }
                else {
                    set temp += [0.];
                }
            }
            set o += [temp];
        }
        return o;
    }
    

    operation U_f(A : Double[][], t : Int, qubits : Qubit[]) : Unit is Ctl {
        mutable eiAt = [[], size = 0];
        for row in A {
            mutable payload = [Complex(0.0, 0.0), size = 0];
            for i in row {
                set payload += [PowC(Complex(E(), 0.0), Complex(0.0, IntAsDouble(t) * i))];
            }
            set eiAt += [payload];
        }
        ApplyUnitary(eiAt, LittleEndian(qubits));
    }


    operation U_f_dag(A : Double[][], t : Int, qubits : Qubit[]) : Unit is Ctl {
        mutable eiAt = [[], size = 0];
        for row in A {
            mutable payload = [Complex(0.0, 0.0), size = 0];
            for i in row {
                set payload += [PowC(Complex(E(), 0.0), Complex(0.0, IntAsDouble(t) * -i))];
            }
            set eiAt += [payload];
        }
        ApplyUnitary(eiAt, LittleEndian(qubits));
    }


    @EntryPoint()
    operation MainOp() : Unit {

        let widthA = 4;  //Should be a power of 2. It is one more than the polynomial degree.
        let data = [[0., 1.], [2., 4.], [3., 5.], [4., 10.], [4., 4.], [7., 3.], [7., 2.], [5., 1.]];

        use b = Qubit[Ceiling(Lg(IntAsDouble(Length(data))))];
        prepareStateB(data, b);    
        DumpMachine();
        


        let Aoriginal = prepareOriginalMatrixA(data, widthA);
        let A = convertAtoHermitian(Aoriginal);
        
        displayMatrix(Aoriginal, "Original A");
        displayMatrix(A, "Hermitian A");




        ResetAll(b);
    }
}
