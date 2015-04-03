# ===============================
# Written by AAE
# TU Delft, Winter 2014
# simulkade.com
# Last edited: 29 December, 2014
# ===============================



function fluxLimiter(flName::ASCIIString)
# This function returns a function handle to a flux limiter of user's
# choice.
# available flux limiters are: 'CHARM', 'HCUS', 'HQUICK', 'VanLeer',
# 'VanAlbada1', 'VanAlbada2', 'MinMod', 'SUPERBEE', 'Sweby', 'Osher',
# 'Koren', 'smart', 'MUSCL', 'QUICK', 'MC', and 'UMIST'. 
# Default limiter is 'SUPERBEE'. See:
# <http://en.wikipedia.org/wiki/Flux_limiter>

# find the flux limiter function

if flName=="CHARM"
  FL(r)=((r>0.0).*r.*(3.0*r+1.0)./(((r+1.0).^2.0)+eps()*(r==-1.0)))
elseif flName=="HCUS"
  FL(r)=(1.5*(r+abs(r))./(r+2.0))
elseif flName=="HQUICK"
  FL(r)=(2.0*(r+abs(r))./((r+3.0)+eps()*(r==-3.0)))
elseif flName=="ospre"
  FL(r)=((1.5*r.*(r+1.0))./(r.*(r+1.0)+1.0+eps()*((r.*(r+1.0)+1.0)==0.0)))
elseif flName=="VanLeer"
  FL(r)=((r+abs(r))./(1.0+abs(r)))
elseif flName=="VanAlbada1"
  FL(r)=((r+r.*r)./(1.0+r.*r))
elseif flName=="VanAlbada2"
  FL(r)=(2.0*r./(1+r.*r))
elseif flName=="MinMod"
  FL(r)=((r>0.0).*min(r,1.0))
elseif flName=="SUPERBEE"
  FL(r)=(max(0.0, max(min(2.0*r,1.0), min(r,2.0))))
elseif flName=="Osher"
  b=1.5
  FL(r)=(max(0.0, min(r,b)))
elseif flName=="Sweby"
  b=1.5
  FL(r)=(max(0.0, max(min(b*r,1.0), min(r,b))))
elseif flName=="smart"
  FL(r)=(max(0.0, min(4.0,min(0.25+0.75*r, 2.0*r))))
elseif flName=="Koren"
  FL(r)=(max(0.0, min(2.0*r, min((1.0+2.0*r)/3.0, 2.0))))
elseif flName=="MUSCL"
  FL(r)=(max(0.0, min(2.0*r, min(0.5*(1+r), 2.0))))
elseif flName=="QUICK"
  FL(r)=(max(0.0, min(2.0, min(2.0*r, (3.0+r)/4.0))))
elseif flName=="UMIST"
  FL(r)=(max(0.0, min(2.0, min(2.0*r, min((1.0+3.0*r)/4.0, (3.0+r)/4.0)))))
else
  println("The flux limiter of your choice is not available. The SUPERBEE flux limiter is used instead.")
  FL(r)=(max(0.0, max(min(2.0*r,1.0), min(r,2.0))))
end

end

function faceEval(f::Function, x::FaceValue)
FaceValue(x.domain,
    f(x.xvalue),
    f(x.yvalue),
    f(x.zvalue))
end

function cellEval(f::Function, x::CellValue)
CellValue(x.domain,
    f(x.value))
end
