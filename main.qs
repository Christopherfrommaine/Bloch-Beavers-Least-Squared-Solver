namespace Least.Squares.Solver {

    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;

    //A little update to test requests

    //TODO:
    //Encode Vectors/Arrays as Amplitudes (For A, b, and similar (A~. and B~))
    //Make the HHL operation
    //Make the QPE/QAE operation
    //Use HHL and QAE to make QsolX
    //Use HHL and QAE and QsolX to make QsolR
    //Use QAE to make norm of X
    //Use QAE to make norm of R
    //Make the f(x) = log(norm(R))^2 + log(norm(X))^2 operation
    //Make the QMF algorithm

    //Specifics of every operation we need to implement:
    
    //Encode Vectors/Arrays as Amplitudes (Aditya)
        //Inputs: datapoints (a classical set of (x, y) points), mu (a binary qubit register representing the mu value), width (a classical value to tell you the degree of the polynomial)
        //Outputs: a qubit register with A~ encoded into the amplitudes, a qubit register with b~ encoded into the amplitudes
        //Example:
        //Given datapoints (0, 1), (2, 3), (5, 3), (7, 1), a mu of 1, and a width of 4
        //b is the y values: [1, 3, 3, 1]
        //A is the x values to the degree of the polynomial:
            //|0^3, 0^2, 0^1, 0^0|   |0,   0,  0, 1|
            //|2^3, 2^2, 2^1, 2^0|   |8,   4,  2, 1|
            //|5^3, 5^2, 5^1, 5^0| = |125, 25, 5, 1|
            //|7^3, 7^2, 7^1, 7^0|   |343, 49, 7, 1|
        
        //Append 0s to the bottom of b to get b~:
            //[1, 3, 3, 1, 0, 0, 0, 0]
        
        //Append mu * I (the identity matrix) to the bottom of A in order to get A~:
            //|0,   0,  0, 1|
            //|8,   4,  2, 1|
            //|125, 25, 5, 1|
            //|343, 49, 7, 1|
            //|1,   0,  0, 0|
            //|0,   1,  0, 0|
            //|0,   0,  1, 0|
            //|0,   0,  0, 1|
        
        //Next, normalize the vector b~ (size of 1) and encode it as the magnitudes of a (in this case) 3-qubit register:
        //1/whatever * (1|000> + 3|001> + 3|010> + 1|011> + 0|100> + 0|101> + 0|110> + 0|111>)

        //Next, do the same for every column of A, and store it as a list of column vectors (Qubit[column index][row index]). (NOTE: it may actually be stores as a list of row vectors, so Qubit[row index][column index]. I don't really know, but be prepared to change it if that is the case.)

        //Return a tuple of (b~, A~)
    

    //Make the HHL Operation (Tony)
        //I'm still a bit unsure of the specific details

    
    //Make the QPE_G Operation (Christopher)
        //Inputs: Qubit register, number of control qubits to use (denoted l in the paper)
        //Output: 
        //Finds the phase of the first element of the register (finds θ given φ).
        //NOTE: The built-in operation "unitary" is needed and extremely helpful here

        //Find a matrix U such that U(|0>) = register (the register you are given)
        //Fing G based on U (formula for G is in paper)
        //Create some control qubits |+>
        //Repeatedly apply powers of G for each control qubit (see paper)
        //Apply the QFT  to the control qubits
        //What you get out is most likely the binary representation of θ*2^(number of control qubits)
        //Return the control qubits
    
    //Make the QAE Operation
        //Inputs: Qubit register, number of control qubits to use (denoted l in the paper)
        //Output: The state |f(cos(theta))>

        //Given a function f (from who knows where. Try to figure that out in the paper)
        //Create U_f (probably with the unitary built in)
        //Apply QPE_G, Apply U_f(control bits, target), Apply adjoint QPE_G
        //Return target
    

    //Make the QMF algorithm
        //See paper (pg 6).
        //Given a mu,




    

        

            







    @EntryPoint()
    operation MainOp() : Result[] {

        use q = Qubit[4];
        ApplyToEach(H, q);
        return MultiM(q);


    }
}
