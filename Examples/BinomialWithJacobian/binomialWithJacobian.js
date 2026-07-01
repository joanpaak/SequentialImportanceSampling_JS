/**
 * y ~ Bernoulli(p)
 * plogis(p) ~ Beta(3, 3)
 */

// Generate data

let n = 300;
let genTheta = -1;
let y = new Matrix({ncol : 1, nrow : n});

for(let i = 0; i < n; i++){
    y.setRow(
        i,
        Stats.rbernoulli(1, Stats.plogis(genTheta, 0, 1))
    )
}

//n = 3;
//let y = new Matrix({ncol : 1, nrow : 3});
//y.setCol(0, [1, 1, 1]);

// Run sequential importance sampling algorithm

let sis = new SIS({
        drawFromPrior : function(n){
            let theta = new Matrix({ncol : 1, nrow : n});

            for(let i = 0; i < n; i++){
                theta.setRow(
                    i, 
                    Stats.qlogis(Stats.rbeta(1, 3, 3)[0], 0, 1));
            }

            return theta;
        },

        prior : function(theta){
            let d = new Array(theta.nrow);

            for(let i = 0; i < theta.nrow; i++){
                d[i] = Stats.dbeta(
                    Stats.plogis(theta.getRow(i), 0, 1), 3, 3) * 
                    Stats.dlogis(theta.getRow(i), 0, 1);
            }

            return d;
        },

        likelihood : function(y, theta){
            let lh = new Array(theta.nrow);

            for(let i = 0; i < theta.nrow; i++){
                lh[i] = Stats.dbernoulli(
                    y, Stats.plogis(theta.getRow(i), 0, 1));
            }

            return lh;
        },

        nParticles : 1000
    },
    {
        logging : true
    }
);

for(let i = 0; i < n; i++){
    sis.addObservation(y.getRow(i));
}

console.log("Generating theta: " + genTheta);
console.log(`Posterior mean after ${n} observations: ${sis.getMarginalMus()}`);
console.log(`Posterior sd after ${n} observations: ${sis.getMarginalSDs()}`);

