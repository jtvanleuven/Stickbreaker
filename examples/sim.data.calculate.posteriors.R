fit.results <- fit.models(Chou.data, c(0.1, 10), 1.1, c(2,1))
sim.data.calculate.posteriors(fit.results$fit.smry
                              ,Chou.data,
                              "Test",
                              c(0.05, 0.5),
                              c(0, 0.25),
                              c(0.1, 10),
                              1.1,
                              c(2,1),
                              50,
                              -1,25)
