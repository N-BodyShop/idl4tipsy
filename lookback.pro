function lbfunc, x
;THIS FUNCTION DEFINES THE INTEGRAND IN THE INTEGRAL THAT HAS TO BE COMPUTED
;IN ORDER TO COMPUTE THE LOOKBACK TIME.  SEE 'THE COSMOLOGICAL CONSTANT' BY
;SEAN CARROLL, EQUATION 36.

omega_m = 0.3
omega_v = 0.7
H = 70.0 * 1.02528 * 10.0^(-12.0) ;in years^-1

return, ((omega_m/x + omega_v * x * x)^(-0.5))/H

end

function lookback, z
;COMPUTES THE LOOKBACK TIME IN YEARS TO AN OBJECT AT REDSHIFT Z

return, qsimp('lbfunc', 1.0/(1.0 + z), 1.0)

end
