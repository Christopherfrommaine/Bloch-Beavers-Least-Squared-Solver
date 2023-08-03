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
    
    

    function createVectorB (equations : Double[][]) : Double[] {
        mutable sumOfSquares = 0.;
        for i in equations {set sumOfSquares += i[Length(equations[0]) - 1];} //Find largest element of b
        let normFactor = Sqrt(sumOfSquares);

        mutable o = [0., size=0];
        for i in equations {
            set o += [i[Length(equations[0]) - 1] / normFactor];
        }
        return o;
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

    
    function prepareMatrixA (equations : Double[][]) : Double[][] {
        mutable o = [[0.], size=0];
        for i in 0 .. Length(equations) - 1 {
            mutable temp = [0., size=0];
            for j in 0 .. Length(equations[i]) - 2 {
                set temp += [equations[i][j]];
            }
            set o += [temp];
        }
        return o;
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
        //Message($"Matrix Form of U: {sum} \nInputs | t: {t}, A: {A}");
        ApplyUnitary(sum, LittleEndian(qubits));
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


        //Configuration
        let totalIterations = 10000;
        let autoGenerateMatrices = false;
        let autoGenerateEigenvalues = false;
        let autoGenerateLoggingAmount = true; //0 is Only Output, 1 is Progress, 2 is Messages, 3 is Final DumpMachine, 4 is Everywhere DumpMachine

        //Inputs
        let linearEquations = [[1., -1./3., 0.], [-1./3., 1., 1.]]; //In the form [a, b, c] such that ax + by = c

        //Direct Input
        let directInput_l = 1;
        let directInput_A = [[1., -1. / 3.], [-1. / 3., 1.]];
        let directInput_bInput = [0., 1.];
        let directInput_cLength = 2;
        let directInput_t = 3. * PI() / 4.;
        let directInput_C = 1.;


        // Program -------
        let A = autoGenerateMatrices ? prepareMatrixA(linearEquations) | directInput_A;
        mutable l = 2;
        mutable t = 0.;
        mutable C = 0.;
        mutable cLength = 0;
        if autoGenerateEigenvalues {
            fail "Not Implemented.";
        }
        else {
            set t = directInput_t;
            set C = directInput_C;
            set cLength = directInput_cLength;
        }
        if autoGenerateLoggingAmount {
            if totalIterations == 1 {
                if autoGenerateEigenvalues {
                    set l = 2;
                }
                else {
                    set l = 3;
                }
            }
            elif totalIterations <= 100 {set l = 2;}
            else {set l = 1;}
        }
        else {
            set l = directInput_l;
        }
        
        if l >= 2 {displayMatrixD(A, "A");}
        if l >= 2 {displayMatrixD([createVectorB(linearEquations)], "b");}

        //HHL Algorithm
        mutable numIterationsTotal = 0;
        mutable outcomes = [];
        mutable previousCompletionPercent = -1;
        repeat {
            //State Prep
            use ancilla = Qubit();
            use c = Qubit[cLength];
            use b = Qubit[Ceiling(Lg(IntAsDouble(Length(linearEquations))))];
            
            if l >= 2 {Message("\nOriginal b (|0>)");}
            if l >= 4 {DumpRegister((), b);}
            PrepareArbitraryStateD(autoGenerateMatrices ? createVectorB(linearEquations) | directInput_bInput, LittleEndian(b));
            if l >= 2 {Message($"\nPrepare b (|1>): {createVectorB(linearEquations)}");}
            if l >= 4 {DumpRegister((), b);}

            //QPE
            //QPE State Prep
            ApplyToEach(H, c);

            //QPE Application of U
            for i in Length(c) - 1 .. -1 .. 0 {
                Controlled U([c[i]], (A, t * IntAsDouble(2 ^ (i)), b));
            }

            if l >= 2 {Message("After QPE");}
            if l >= 4 {DumpMachine();}

            Adjoint QFTLE(LittleEndian(c));

            if l >= 2 {Message("After QFT");}
            if l >= 4 {DumpMachine();}

            //QPE Controlled Rotation
            ancillaRotations(c, ancilla);

            if l >= 2 {Message("After Rotation");}
            if l >= 4 {DumpMachine();}

            //Continuing the Algorithm
            //QPE IQPE
            QFTLE(LittleEndian(c));
            for i in Length(c) - 1 .. -1 .. 0 {
                Controlled U([c[i]], (A, -1. * t * IntAsDouble(2 ^ (i)), b));
            }
            ApplyToEach(H, c);

            if l >= 2 {Message("Final State:");}
            if l >= 3 {DumpMachine();}

            if M(ancilla) == Zero {
                set outcomes += [ResultArrayAsInt(MultiM(b))];
                set numIterationsTotal += 1;
            }
            ResetAll([ancilla] + b + c);
            if l >= 1 {
                if numIterationsTotal % Ceiling(IntAsDouble(totalIterations) / 100.) == 0 {
                    mutable progress = "";
                    mutable spaces = "";
                    let percentComplete = (100 * numIterationsTotal) / (totalIterations);
                    if percentComplete <= 9 {set spaces += " ";}
                    if percentComplete <= 99 {set spaces += " ";}
                    for j in 0 .. 20 {
                        set progress +=  100 * j / 20 <= (100 * numIterationsTotal) / (totalIterations) ? "#" | "-";
                    }
                    if percentComplete != previousCompletionPercent {
                        Message($"{(100 * numIterationsTotal) / (totalIterations)}%{spaces} | {progress} |");
                    }
                    set previousCompletionPercent = percentComplete;
                }
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



        if l >= 2 {Message($"Outcomes: {outcomes}");}
        if outcomesFinal[0] < outcomesFinal[1] {
            Message($"\n--------\nOutput: {outcomesFinal} | 1 : {IntAsDouble(outcomesFinal[1]) / IntAsDouble(outcomesFinal[0])}\n--------\n");
        }
        else {
            Message($"\n--------\nOutput: {outcomesFinal} | {IntAsDouble(outcomesFinal[0]) / IntAsDouble(outcomesFinal[1])} : 1\n--------\n");
        }
    }
}
