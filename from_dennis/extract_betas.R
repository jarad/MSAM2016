library(dplyr)

d = readr::read_csv("BdrawsIM.csv")[101:10100,-1]

hartley    = d[,11+c(0,41)]
names(hartley)    = c("intercept","distance")
hartley$fkicker = "GH-060"

janikowski = d[,40+c(0,41)]
names(janikowski) = c("intercept","distance")
janikowski$fkicker = "SJ-030"

d = bind_rows(hartley, janikowski)

predict = function(d) {
  expand.grid(dist = 17:80,
             synth = c(0,1)) %>%
    mutate(probability = pnorm(d$intercept + d$distance*dist + 0.105*synth))
}

p_d = d %>%
  tidyr::gather(parameter, draw, intercept, distance) %>%
  group_by(fkicker, parameter) %>%
  summarize(mean = mean(draw)) %>%
  tidyr::spread(parameter, mean) %>%
  do(predict(.))

ggplot(p_d, aes(dist, probability, color=synth, group=synth)) +
  geom_line() + 
  facet_grid(fkicker~.)