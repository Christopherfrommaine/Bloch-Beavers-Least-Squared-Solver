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



    //Prep (Almost Done)
        //Prep b
        //Prep A (Hermitian)
    //Find Eigenvectors of A
    //Get U given A an t
        //U = e^iAt
        //Convert A to Adiag
        //Find Udiag
        //Convert to U

        //V = (u0, u1) (The Eigenvectors)
        //V_t = daggar of V = V
        //Adiag = V_t A V
        //Udiag = (e^i eigenvalue0 t   0                )
        //        (0                   e^i eigenvalue1 t)
        //U = V Udiag V_t
        //U^2 = V (Udiag^2) V_t
        //Etc for all powers

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
    
    

    
    

    operation prepareStateB (data : Double[][], qubits : Qubit[], entangledAmplitudeb : Qubit) : Unit {
        //This prepares the vector |b> as the amplitudes of the |1> state of an ancillary qubit
        //when entangled with the index qubit. ie: b[5] = |5>0.7|1>
        //This is the method used in the 25 page paper.

        mutable maxB = 0.;
        for i in data {if i[1] > maxB {set maxB = i[1];}} //Find largest element of b

        mutable b = [];
        for i in data {set b += [i[1] / maxB];} //Gets the b vecotor normalized w/ all entries < 1

        let n_b = Ceiling(Lg(IntAsDouble(Length(b)))); //Find qubit length of b

        ApplyToEach(H, qubits);
        for i in 0 .. 2 ^ n_b - 1 {
            let binaryrepresentation = IntAsBoolArray(i, n_b);
            ApplyPauliFromBitString(PauliX, false, binaryrepresentation, qubits);
            Controlled Ry(qubits, (2. * ArcSin(b[i]), entangledAmplitudeb));
            ApplyPauliFromBitString(PauliX, false, binaryrepresentation, qubits);
        }
    }


    function displayMatrix(matrix : Double[][]) : Unit {
        mutable o = "|";
        for column in matrix {
            for element in column {
                set o += DoubleAsString(element) + " ";
            }
            set o += "|\n|";
        }
        Message($"Matrix: {o}");
    }


    function transpose(matrix : Double[][]) : Double[][] {
        mutable o = [[0.], size=0];
        for i in 0 .. Length(matrix[0]) - 1 {
            mutable temp = [0., size=0];
            for j in 0 .. Length(matrix) - 1 {
                set temp += [matrix[i][j]];
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
                if i > h and j <= h {
                    set temp += [A[i - h][j]];
                }
                elif i <= h and j > h {
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





    @EntryPoint()
    operation MainOp() : Unit {

        let widthA = 4;  //Should be a power of 2. It is one more than the polynomial degree.
        let data = [[0., 1.], [2., 4.], [3., 5.], [4., 10.], [4., 4.], [7., 3.], [7., 2.], [5., 1.]];


        use (bAmplitude, bIndex) = (Qubit(), Qubit[Ceiling(Lg(IntAsDouble(Length(data))))]);
        prepareStateB(data, bIndex, bAmplitude);
        
        DumpMachine();
        ResetAll(bIndex + [bAmplitude]);


        let Aoriginal = prepareOriginalMatrixA(data, widthA);
        Message($"Aoriginal: {Aoriginal}");
        let A = convertAtoHermitian(Aoriginal);
        Message($"Ahermitian: {A}");
        





    }
}
