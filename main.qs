namespace Least.Squares.Solver {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;


    // @EntryPoint() denotes the start of program execution.
    @EntryPoint()
    operation MainOp() : Int {
        
        use register = Qubit[8];

        ApplyToEach(H, register);

        return ResultArrayAsInt(MultiM(register));

    }
}
