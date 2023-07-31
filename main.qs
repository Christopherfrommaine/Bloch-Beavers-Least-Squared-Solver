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
        mutable o = name;

        mutable maxLength = 0;
        for column in matrix {
            for element in column {
                if element > 1. and (Ceiling(Log10(element))) > maxLength {set maxLength = Ceiling(Log10(element));}
            }
        }

        for column in matrix {
            set o += "\n|";
            for element in column {
                set o += DoubleAsString(element);
                for space in 0 .. maxLength - (element > 1. ? Ceiling(Log10(element)) | 1) {set o += " ";}
            }
            set o += "|";
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
    

    operation U_f(A : Double[][], t : Double, qubits : Qubit[]) : Unit is Ctl {
        mutable eiAt = [[], size = 0];
        for row in A {
            mutable payload = [Complex(0.0, 0.0), size = 0];
            for i in row {
                set payload += [PowC(Complex(E(), 0.0), Complex(0.0, t * i))];
            }
            set eiAt += [payload];
        }
        ApplyUnitary(eiAt, LittleEndian(qubits));
    }


    operation U_f_inverse(A : Double[][], t : Double, qubits : Qubit[]) : Unit is Ctl {
        mutable eiAt = [[], size = 0];
        for row in A {
            mutable payload = [Complex(0.0, 0.0), size = 0];
            for i in row {
                set payload += [PowC(Complex(E(), 0.0), Complex(0.0, t * -i))];
            }
            set eiAt += [payload];
        }
        ApplyUnitary(eiAt, LittleEndian(qubits));
    }


    operation QCR(b : Qubit[], c : Qubit[], U : (Qubit[] => Unit is Ctl)) : Unit{
        for controlQubit in c {
            Controlled U([controlQubit], b);
        }
    }

    operation ancillaRotations(c : Qubit[], a : Qubit) : Unit {
        for i in 0 .. Length(c) - 1 {
            Controlled Ry([c[i]], (2. * ArcCos(1. / (IntAsDouble(i) + 1.)), a));
        }
    }


    @EntryPoint()
    operation MainOp() : Unit {

        //Major Steps:
            //State Prep
            //QPE
                //State Prep (ApplyToEach(H, Controls))
                //Application of U
                //IQFT
            //Controlled Ry onto Ancilla
            //Measure Ancilla (If 1, continue. If 0, restart)
            //IQPE
            //Measure b register to find x


        //Inputs
        
        let inputMatrixFormsDirectly = true;
        

        let widthA = 8;  //Only specific numbers work. It is a little weird. Ask me (Christopher) to explain it if this causes problems later
        let data = [[0., 1.], [2., 4.], [3., 5.], [4., 10.], [4., 4.], [7., 3.], [7., 2.], [5., 1.]];

        let Aoriginal = prepareOriginalMatrixA(data, widthA);
        

        //State Preperation (Need to define variables outside of scope of repeat loop)
        if not inputMatrixFormsDirectly {displayMatrix(Aoriginal, "Original A");}
        

        let directInput_A = [[1., -1. / 3.], [-1. / 3., 1.]];
        let directInput_bInput = [0., 1.];
        let directInput_bLength = 1;

        let A = inputMatrixFormsDirectly ? directInput_A | convertAtoHermitian(Aoriginal);
        displayMatrix(A, "A");

        use b = Qubit[inputMatrixFormsDirectly ? directInput_bLength | Ceiling(Lg(IntAsDouble(Length(data)))) + 1];

        if not inputMatrixFormsDirectly {
            prepareStateB(data, b);
        }
        else {
            PrepareArbitraryStateD(directInput_bInput, LittleEndian(b));
            DumpMachine();
        }
        
        

        //(For QPE)
        use c = Qubit[10];
        use ancilla = Qubit();

        let t = 3. * PI() / 4.;

        //HHL Algorithm
        mutable dontRepeatComputation = false;
        mutable repitionCount = 0;
        repeat {
            //State Prep Continued
            if not inputMatrixFormsDirectly {
                prepareStateB(data, b);
            }
            else {
                PrepareArbitraryStateD(directInput_bInput, LittleEndian(b));
            }   
            //DumpMachine();


            //QPE
            //QPE State Prep
            ApplyToEach(H, c);

            //QPE Application of U
            QCR(b, c, U_f(A, t, _)); //Im slightly unsure of the naming if this should actually be called QCR. It should work though

            Adjoint QFT(LittleEndianAsBigEndian(LittleEndian(c)));

            //QPE Controlled Rotation
            ancillaRotations(c, ancilla);

            set dontRepeatComputation = M(ancilla) == One;
            if not dontRepeatComputation {
                ResetAll(b + c);
            }
            Message($"Ancilla measured to be " + (dontRepeatComputation ? "One. Continuing..." | "Zero Repeating..."));
            set repitionCount += 1;
        }
        until dontRepeatComputation or repitionCount > 100;
        if repitionCount > 100 {Message("Computation Failed: Ancilla never measured to be One");}

        //QPE IQPE
        QFT(LittleEndianAsBigEndian(LittleEndian(c)));
        QCR(b, c, U_f_inverse(A, t, _));
        ApplyToEach(H, c);

        //Reset
        ResetAll(c);

        //Final Measureement
        DumpRegister((), b);
        Message($"Output: {MultiM(b)}");
    }
}
