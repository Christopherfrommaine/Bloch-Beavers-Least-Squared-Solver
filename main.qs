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

    

    operation prepareStateB (data : Int[][], qubits : Qubit[], entangledAmplitudeb : Qubit) : Unit {
        //This prepares the vector |b> as the amplitudes of the |1> state of an ancillary qubit
        //when entangled with the index qubit. ie: b[5] = |5>0.7|1>
        //This is the method used in the 25 page paper.

        mutable maxB = 0.;
        for i in data {if IntAsDouble(i[1]) > maxB {set maxB = IntAsDouble(i[1]);}} //Find largest element of b

        mutable b = [];
        for i in data {set b += [IntAsDouble(i[1]) / maxB];} //Gets the b vecotor normalized w/ all entries < 1

        let n_b = Ceiling(Lg(IntAsDouble(Length(b)))); //Find qubit length of b

        ApplyToEach(H, qubits);
        for i in 0 .. 2 ^ n_b - 1 {
            let binaryrepresentation = IntAsBoolArray(i, n_b);
            ApplyPauliFromBitString(PauliX, false, binaryrepresentation, qubits);
            Controlled Ry(qubits, (2. * ArcSin(b[i]), entangledAmplitudeb));
            ApplyPauliFromBitString(PauliX, false, binaryrepresentation, qubits);
        }
    }





    @EntryPoint()
    operation MainOp() : Unit {

        let polynomialDegree = 3;  //Should be a power of 2 minus 1
        let data = [[0, 1], [2, 4], [3, 5], [4, 10], [4, 4], [7, 3], [7, 2], [5, 1]];


        use (bAmplitude, bIndex) = (Qubit(), Qubit[Ceiling(Lg(IntAsDouble(Length(data))))]);
        prepareStateB(data, bIndex, bAmplitude);
        
        DumpMachine();
        ResetAll(bIndex + [bAmplitude]);


        

    }
}
