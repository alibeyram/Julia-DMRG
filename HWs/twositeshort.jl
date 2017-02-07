sz = Float64[0.5 0; 
	     0 -0.5]
sp = Float64[0 1; 
	     0 0]
sm = sp'
one = eye(2)

Htwosite = 0.5 * (kron(sp,sm) + kron(sm,sp)) + kron(sz,sz)
@show size(Htwosite)
Hmat = reshape(Htwosite,4,4)
ev2 = eigfact(Hmat)
@show ev2[:values]

