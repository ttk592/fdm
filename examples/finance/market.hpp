#ifndef MF_FDM_MARKET_HPP
#define MF_FDM_MARKET_HPP

namespace fin
{

struct MarketHeston {
    double	spot, vol;	// initial spot and vol
    double	rd, rf;		// domestic, foreign interest rate
    double	theta;		// longterm variance average
    double	kappa;		// strength of mean reversion
    double	xi;		// vol of vol
    double	rho;		// correlation between spot and vol
};

} // namespace fin

#endif /* MF_FDM_MARKET_HPP */
