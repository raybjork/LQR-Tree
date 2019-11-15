function Sdot = TV_Ricatti(t, S, system)
    S = reshape(S, [2, 2]);
    A = system.A;
    B = system.B;
    Q = system.Q;
    R = system.R;

    Sdot = S*B*inv(R)*B'*S - Q - A'*S - S*A;
    Sdot = reshape(Sdot, [4,1]);
end