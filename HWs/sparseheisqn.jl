const n = 20
const sztot = 0			# 2 * Sztot, so an integer
@show (n,sztot)

const si1 = (fill(2,n)...)	# array of size n, filled with 2's, convert to tuple

function mysub2ind(inds::Vector{Int64},n::Int64)	# This tests much faster than sub2ind
    res = inds[n]-1
    for i=n-1:-1:1
	res *= 2
	res += inds[i]-1
    end
    res+1
end

function getjind(inds::Vector{Int64},j)	# get the index based on inds, but with j and j+1 sites swapped
    (inds[j],inds[j+1]) = (inds[j+1],inds[j])
    jind = mysub2ind(inds,n)
    (inds[j],inds[j+1]) = (inds[j+1],inds[j])	# put back in order, since array not copied
    jind
end

function makeblock()
    inblock = Vector{Int64}(0)
    for i=1:2^n
	inds = [ind2sub(si1,i)...]
	sz = sum(j -> 2*inds[j]-3,1:n)
	sz == sztot && push!(inblock,i)
    end
    inblock
end

@time inblock = makeblock()
const nb = length(inblock)
@show nb

@time tuplist = [(inblock[i],i) for i=1:nb]
@time maptoind = Dict{Int64,Int64}(tuplist)
println("Done setting up map")
flush(STDOUT)

# Constructing the sparse matrix from three lists (row,col,value) is MUCH faster
# than using spzeros to make an empty sparse matrix, and assigning to elements.
# for n=20, about 100 times faster!
    
row = Vector{Int64}(0)
col = Vector{Int64}(0)
val = Vector{Float64}(0)

diagterm(a) = 0.25 * sum(j->(a[j] == a[j+1] ? 1 : -1) ,1:length(a)-1)

function dopush!(i,j,v) 	# put this element in terms lists
    push!(row,i)
    push!(col,j)
    push!(val,v)
end
function getHlist()
    for i=1:nb
	ib = inblock[i]
	inds = [ind2sub(si1,ib)...]
	dopush!(i,i,diagterm(inds))
	for j=1:n-1
	    inds[j] == inds[j+1] && continue
	    jind=getjind(inds,j)
	    dopush!(i,maptoind[jind],0.5)
	end
    end
end
@time getHlist()
@show length(row)
println("Done making H list")
flush(STDOUT)

@time H = sparse(row,col,val)
println("Done making H")
flush(STDOUT)

@time evn = eigs(H;nev=1, which=:SR)
println("Done with Lanczos")
flush(STDOUT)

@show evn[1]
@show evn[3:5]
