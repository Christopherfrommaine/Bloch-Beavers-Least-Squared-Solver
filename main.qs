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


    function displayMatrixD(inputMatrix : Double[][], name : String) : Unit {
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


    function matMulD(m1 : Double[][], m2 : Double[][]) : Double[][] {
        if Length(m1) != Length(m2[0]) {Message("ERROR! Incompatable Matrices in matMul");}

        //ONLY WORKS ON SQUARE MATRICES

        mutable o = [[0.], size=0];
        for column in 0 .. Length(m1) - 1 {
            mutable columnTemp = [0., size=0];
            for row in 0 .. Length(m1[column]) - 1 {
                mutable sum = 0.;
                for i in 0 .. Length(m1) - 1 {
                    set sum += (m1[i][row]) * (m2[column][i]);
                }
                set columnTemp += [sum];
            }
            set o += [columnTemp];
        }
        return o;
    }
    function matMulC(m1 : Complex[][], m2 : Complex[][]) : Complex[][] {
        if Length(m1) != Length(m2[0]) {Message("ERROR! Incompatable Matrices in matMul");}

        //ONLY WORKS ON SQUARE MATRICES

        mutable o = [[Complex(0., 0.)], size=0];
        for column in 0 .. Length(m1) - 1 {
            mutable columnTemp = [Complex(0., 0.), size=0];
            for row in 0 .. Length(m1[column]) - 1 {
                mutable sum = Complex(0., 0.);
                for i in 0 .. Length(m1) - 1 {
                    set sum = PlusC(sum, TimesC(m1[i][row], m2[column][i]));
                }
                set columnTemp += [sum];
            }
            set o += [columnTemp];
        }
        return o;
    }


    function scalarMatMulC(m1 : Complex[][], scalar : Complex) : Complex[][] {
        mutable o = [[Complex(0., 0.)], size=0];
        for i in 0 .. Length(m1) - 1 {
            mutable temp = [Complex(0., 0.), size=0];
            for j in 0 .. Length(m1[i]) - 1 {
                set temp += [TimesC(scalar, m1[i][j])];
            }
            set o += [temp];
        }
        return o;
    }


    function AddMatC(m1 : Complex[][], m2 : Complex[][]) : Complex[][] {
        mutable o = [[Complex(0., 0.)], size=0];
        for i in 0 .. Length(m1) - 1 {
            mutable temp = [Complex(0., 0.), size=0];
            for j in 0 .. Length(m1[i]) - 1 {
                set temp += [PlusC(m1[i][j], m2[i][j])];
            }
            set o += [temp];
        }
        return o;
    }


    function isMatrixEqual (m1 : Double[][], m2 : Double[][]) : Bool {
        if Length(m1) != Length(m2) or Length(m1[0]) != Length(m2[0]) {return false;}
        for i in 0 .. Length(m1) - 1 {
            for j in 0 .. Length(m1[i]) - 1 {
                if m1[i][j] != m2[i][j] {
                    return false;
                }
            }
        }
        return true;
    }

    function doubleMatrixToComplex(m : Double[][]) : Complex[][] {
        mutable o = [[Complex(0., 0.)], size=0];
        for i in 0 .. Length(m) - 1 {
            mutable temp = [Complex(0., 0.), size=0];
            for j in 0 .. Length(m[i]) - 1 {
                set temp += [Complex(m[i][j], 0.)];
            }
            set o += [temp];
        }
        return o;
    }

    function IdentityI(size : Int) : Int[][] {
        mutable o = [[0], size=0];
        for i in 0 .. size - 1 {
            mutable temp = [0, size=0];
            for j in 0 .. size - 1 {
                set temp += [i == j ? 1 | 0];
            }
            set o += [temp];
        }
        return o;
    }
    function IdentityD(size : Int) : Double[][] {
        mutable o = [[0.], size=0];
        for i in 0 .. size - 1 {
            mutable temp = [0., size=0];
            for j in 0 .. size - 1 {
                set temp += [i == j ? 1. | 0.];
            }
            set o += [temp];
        }
        return o;
    }
    function IdentityC(size : Int) : Complex[][] {
        mutable o = [[Complex(0., 0.)], size=0];
        for i in 0 .. size - 1 {
            mutable temp = [Complex(0., 0.), size=0];
            for j in 0 .. size - 1 {
                set temp += [i == j ? Complex(1., 0.) | Complex(0., 0.)];
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


    operation U(A : Double[][], tTemp : Double, qubits : Qubit[]) : Unit is Ctl {
        let t = -1. * tTemp;
        //Inverse can be taken by setting t to be negativec
        //Powers can be taken by multiplying t
        
        let iterations = 20;
        mutable sum = [[Complex(0., 0.), Complex(0., 0.)], [Complex(0., 0.), Complex(0., 0.)]];
        for n in 0 .. iterations {
            //let s = PowC(TimesC(Complex(t / IntAsDouble(FactorialI(n)), 0.), Complex(0., 1.)), Complex(IntAsDouble(n), 0.));
            let s = TimesC(PowC(Complex(0., -1. * t), Complex(IntAsDouble(n), 0.)),  Complex(1. / IntAsDouble(FactorialI(n)), 0.)); //Computes (-i * t) ^ n  / n!
            mutable m = IdentityC(Length(A));
            for _ in 1 .. n {
                set m = matMulC(m, doubleMatrixToComplex(A));
            }
            set m = scalarMatMulC(m, s);
            set sum = AddMatC(sum, m);
        }
        if false {
            mutable normFactor = 0.;
            for i in sum {
                for j in i {
                    set normFactor += AbsComplex(j) ^ 2.;
                }
            }
            set normFactor = 1. / Sqrt(normFactor);
            mutable o = [[Complex(0., 0.)], size=0];
            for i in sum {
                mutable temp = [Complex(0., 0.), size=0];
                for j in i {
                    set temp += [TimesC(j, Complex(normFactor, 0.))];
                }
                set o += [temp];
            }
            Message($"Normalized Matrix Form of U: {o} \nInputs | t: {t}, A: {A}");
            ApplyUnitary(o, LittleEndian(qubits));
        }
        else {
            Message($"Matrix Form of U: {sum} \nInputs | t: {t}, A: {A}");
            ApplyUnitary(sum, LittleEndian(qubits));
        }
    }


    operation ancillaRotations(c : Qubit[], a : Qubit) : Unit {
        for i in 0 .. Length(c) - 1 {
            Controlled Ry([c[i]], (2. * ArcCos(1. / (IntAsDouble(i) + 1.)), a));
        }
    }

    

    @EntryPoint()
    operation MainOp() : Unit {
        
        Message("Starting HHL Algorithm \n\n\n\n\n\n\n\n\n\n\n\n\n");

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

        let totalIterations = 1;
        mutable numIterationsTotal = 0;
        mutable outcomes = [];
        
        let widthA = 8;  //Only specific numbers work. It is a little weird. Ask me (Christopher) to explain it if this causes problems later
        let data = [[0., 1.], [2., 4.], [3., 5.], [4., 10.], [4., 4.], [7., 3.], [7., 2.], [5., 1.]];

        //State Preperation (Need to define variables outside of scope of repeat loop)
        let Aoriginal = prepareOriginalMatrixA(data, widthA);
        if not inputMatrixFormsDirectly {displayMatrixD(Aoriginal, "Original A");}
        

        let directInput_A = [[1., -1. / 3.], [-1. / 3., 1.]];
        let directInput_bInput = [0., 1.];
        let directInput_bLength = 1;
        
        let A = inputMatrixFormsDirectly ? directInput_A | convertAtoHermitian(Aoriginal);
        displayMatrixD(A, "A");

        let t = 3. * PI() / 4.;

        //HHL Algorithm
        repeat {
            set numIterationsTotal += 1;
            mutable dontRepeatComputation = false;
            mutable repitionCount = 0;
            repeat {
                
                //State Prep
                use b = Qubit[inputMatrixFormsDirectly ? directInput_bLength | Ceiling(Lg(IntAsDouble(Length(data)))) + 1];
                use c = Qubit[2];
                use ancilla = Qubit();
                Message("\nOriginal b (|0>)");
                DumpRegister((), b);
                if not inputMatrixFormsDirectly {prepareStateB(data, b);}
                else {PrepareArbitraryStateD(directInput_bInput, LittleEndian(b));}
                Message("\nPrepare b (|1>)");
                DumpRegister((), b);

                //QPE
                //QPE State Prep
                ApplyToEach(H, c);

                //QPE Application of U
                for i in 0 .. Length(c) - 1 {
                    Controlled U([c[i]], (A, t * IntAsDouble(2 ^ (i)), b));
                }

                Message("After QPE");
                DumpMachine();

                Adjoint QFT(LittleEndianAsBigEndian(LittleEndian(c)));

                Message("After QFT");
                DumpMachine();

                //QPE Controlled Rotation
                ancillaRotations(c, ancilla);

                Message("After Rotation");
                DumpMachine();

                //Ancilla Measurement
                set dontRepeatComputation = M(ancilla) == One;
                Message($"Ancilla measured to be " + (dontRepeatComputation ? "One. Continuing..." | "Zero Repeating..."));
                set repitionCount += 1;

                if dontRepeatComputation {
                    //Continuing the Algorithm
                    //QPE IQPE
                    QFT(LittleEndianAsBigEndian(LittleEndian(c)));
                    for i in 0 .. Length(c) - 1 {
                        Controlled U([c[i]], (A, -1. * t * IntAsDouble(2 ^ (i)), b));
                    }
                    ApplyToEach(H, c);

                    //Final Measureement
                    DumpRegister((), b);
                    DumpMachine();
                    let output = MultiM(b);
                    Message($"Output: {output}");
                    set outcomes += [ResultArrayAsInt(output)];

                    //Reset
                    ResetAll(c);
                }
                else {
                    ResetAll(b + c);
                }
            }
            until dontRepeatComputation; // or repitionCount > 10;
            if repitionCount > 10 and false {
                Message("Computation Failed: Ancilla never measured to be One");
                fail "because why not";
            }
        }
        until numIterationsTotal == totalIterations;

        //Statistical Post-Processing
        mutable outcomesFinal = [];
        mutable maxIndex = 1;
        for i in outcomes {if i > maxIndex {set maxIndex = i;}}

        for i in 0 .. maxIndex {
            mutable count = 0;
            for j in outcomes {
                if j == i {
                    set count += 1;
                }
            }
            set outcomesFinal += [count];
        }



        Message($"Outcomes: {outcomes}");
        if outcomesFinal[0] < outcomesFinal[1] {
            Message($"\n--------\nOutput: {outcomesFinal} | 1 : {IntAsDouble(outcomesFinal[1]) / IntAsDouble(outcomesFinal[0])}");
        }
        else {
            Message($"\n--------\nOutput: {outcomesFinal} | {IntAsDouble(outcomesFinal[0]) / IntAsDouble(outcomesFinal[1])} : 1");
        }
    }
}
