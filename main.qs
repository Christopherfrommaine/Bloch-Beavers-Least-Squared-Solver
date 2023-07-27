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
    

    
    operation prepareAmplitudeB (b : Double[], register : Qubit[]) : Unit{
        mutable ampArray = [b];
        for i in 0 .. Length(register) - 1 {
            mutable tempArray = [0., size=0];
            mutable tempValue = 0.;
            for j in 0 .. Length(ampArray[i]) - 1 {
                set tempValue += ampArray[i][j] ^ 2.;
                if j % 2 == 1 {
                    set tempArray += [Sqrt(tempValue)];
                    set tempValue = 0.;
                }
            }
            set ampArray += [tempArray];
        }
        Message($"Amp: {ampArray}");
        for i in Length(ampArray) - 2 .. -1 .. 0 {
            for j in 0 .. Length(ampArray[i]) - 1 {
                Message($"{i}, {j}");
            }
        }

        for i in Length(ampArray) - 2 .. -1 .. 0 {
            for j in 0 .. Length(ampArray[i]) - 1 {
                let value = ampArray[i][j];
                let theta = 2. * ArcSin(value);
                let binrepr = IntAsBoolArray(j, Ceiling(Lg(IntAsDouble(Length(ampArray[i])))));
                let controlls = i + 1 .. Length(ampArray[i]) - 1;
                Message($"{i}, {j}, {binrepr}");
                Message($"{Length(register)}");
                Message($"{([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])[controlls]}");

                ApplyPauliFromBitString(PauliX, false, binrepr, register[controlls]);
                Controlled Ry(register[controlls], (theta, register[i]));
                ApplyPauliFromBitString(PauliX, false, binrepr, register[controlls]);
            }
        }

    }
    

    operation prepareStateB (data : Double[][], register : Qubit[]) : Unit {
        //This prepares the vector |b> as the amplitudes of the |1> state of an ancillary qubit
        //when entangled with the index qubit. ie: b[5] = |5>0.7|1>
        //This is the method used in the 25 page paper.

        mutable sumOfSquares = 0.;
        for i in data {set sumOfSquares += i[1] ^ 2.;} //Find largest element of b
        let normFactor = Sqrt(sumOfSquares);

        mutable b = [];
        for i in data {set b += [i[1] / normFactor];} //Gets the b vecotor normalized w/ all entries < 1

        Message($"B (normalized): {b}");

        for i in Length(register) - 1 .. -1 .. 0 {
            mutable squareSum = 0.;
            for j in 0 .. Length(b) - 1 {
                if j % 2 ^ i == 1 {
                    set squareSum += b[j] ^ 2.;
                }
            }
            let theta = 2. * ArcSin(Sqrt(squareSum));
            let binrepr = IntAsBoolArray()
            Message($"{binrepr}");
            Message($"{i}, Theta: {theta}, {([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])[i + 1 .. Length(register) - 1]}");
            Controlled Ry(register[i + 1 .. Length(register) - 1], (theta, register[i]));
        }
    }


    function displayMatrix(inputMatrix : Double[][]) : Unit {
        let matrix = transpose(inputMatrix);
        mutable o = "|";

        mutable maxLength = 0;
        for column in matrix {
            for element in column {
                if element > 1. and (Ceiling(Log10(element))) > maxLength {set maxLength = Ceiling(Log10(element));}
            }
        }
        // Message($"MaxLength: {maxLength}");

        for column in matrix {
            for element in column {
                set o += DoubleAsString(element);
                for space in 0 .. maxLength - (element > 1. ? Ceiling(Log10(element)) | 1) {set o += " ";}
            }
            set o += "|\n|";
        }
        Message($"Matrix: \n{o}");
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



    @EntryPoint()
    operation MainOp() : Unit {

        let widthA = 4;  //Should be a power of 2. It is one more than the polynomial degree.
        let data = [[0., 1.], [2., 4.], [3., 5.], [4., 10.], [4., 4.], [7., 3.], [7., 2.], [5., 1.]];


        use b = Qubit[Ceiling(Lg(IntAsDouble(Length(data))))];
        prepareStateB(data, b);
        
        DumpMachine();
        ResetAll(b);


        let Aoriginal = prepareOriginalMatrixA(data, widthA);
        Message($"Aoriginal: {Aoriginal}");
        let A = convertAtoHermitian(Aoriginal);
        Message($"Ahermitian: {A}");

        
        displayMatrix(Aoriginal);
        displayMatrix(A);
        





    }
}
