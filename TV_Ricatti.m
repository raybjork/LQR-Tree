function Ldot = TV_Ricatti(t, S, system)
    S = reshape(S, [2, 2]);
    A = system.A;
    B = system.B;
    Q = system.Q;
    R = system.R;

    Ldot = S*B*inv(R)*B'*S - Q - A'*S - S*A;
    Ldot = reshape(Sdot, [4,1]);
end