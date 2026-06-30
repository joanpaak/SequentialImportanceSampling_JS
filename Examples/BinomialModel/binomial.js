/**
 * BINOMIAL MODEL
 * 
 * y ~ Bernoulli(p)
 * p ~ Beta(3, 3)
 * 
 */


function drawFromPrior(n) {
    let x = Stats.rbeta(n, 3, 3);

    let theta = new Matrix({ ncol: 1, nrow: n });
    theta.setCol(0, x);

    return theta;
}

function prior(theta){
    let ptheta = new Array(theta.nrow);

    for(let i = 0; i < theta.nrow; i++){
        ptheta[i] = Stats.dbeta(theta.getRow(i), 3, 3); 
    }

    return ptheta;
}

function likelihood(y, theta){
    let lh = new Array(theta.nrow);

    for(let i = 0; i < theta.nrow; i++){
        lh[i] = Stats.dbernoulli(y, theta.getRow(i));
    }

    return lh;
}

let sis = new SIS(
    {
        drawFromPrior : drawFromPrior,
        prior : prior,
        likelihood : likelihood,
        nParticles : 1000
    }
);

let n = 140;
let y = Stats.rbernoulli(n, 0.8);

for(let i = 0; i < n; i++){
    sis.addObservation([y[i]]);
}

console.log("Generating value for parameter p: 0.8");
console.log("Posterior mean for parameter p:");
console.log(sis.getMarginalMus());
