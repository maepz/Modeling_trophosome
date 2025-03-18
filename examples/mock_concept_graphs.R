library(ggplot2)

#### 1 ####
time=seq(0,100)
rich_up=seq(0,100)*0.1
rich_down=rev(rich_up)
rich_steady=rep(5,100)

df1=cbind(time,rich_up,rich_down,rich_steady) %>% as.data.frame() %>% reshape2::melt(id.vars=c('time'),variable.name='hp')
df1 %>% head

ggplot(df1)+
  geom_line(aes(x=time,y=value,color=hp))+
  geom_point(data=df1[df1$time %in% c(0,100),],aes(x=time,y=value,color=hp),size=5)+
  scale_x_continuous(limits=c(-20,120))+
  scale_y_continuous(limits=c(-0.5,10))+
  theme_bw()+
  labs(x='time (host generations)',y='richness (environmental symbionts)')+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggplot(df1)+
  geom_line(aes(x=time,y=value,color=hp))+
  # geom_point(data=df1[df1$time %in% c(0,100),],aes(x=time,y=value,color=hp),size=5)+
  # scale_x_continuous(limits=c(-20,120))+
  scale_y_continuous(limits=c(-1,11))+
  scale_color_discrete(guide='none',)+
  theme_bw()+
  labs(x='time (host generations)',y='richness (environmental symbionts)')+
  theme(
    # axis.text.x = element_blank()
        axis.text.y = element_blank())

#### 2 ####
time=seq(0,100)
rich_up=(seq(0,100)*0.025)+5
# rich_down=rev(rich_up)
rich_steady=rep(5,101)

df1=cbind(time,rich_up,rich_steady) %>% as.data.frame() %>% reshape2::melt(id.vars=c('time'),variable.name='hp')
df1 %>% head

ggplot(df1)+
  geom_line(aes(x=time,y=value,linetype=hp))+
  # geom_point(data=df1[df1$time %in% c(0,100),],aes(x=time,y=value,color=hp),size=5)+
  # scale_x_continuous(limits=c(-20,120))+
  scale_y_continuous(limits=c(2,11))+
  # scale_color_discrete(guide='none',)+
  scale_linetype_discrete(labels=c('with symbiosis','without symbiosis'),name='')+
  theme_bw()+
  labs(x='Time (Generations)',y='Symbiont richness')+
  theme(legend.position="top",
    # axis.text.x = element_blank()
    axis.text.y = element_blank())

#### ####
my_data <- data.frame(mean =  c(0, 0, 0,0,2,-2,0,0),
                      stdev = c(1, 1, 1, 1,1,1,2,0.5), 
                      test = c("H1", "H2", "H3", "H4","H1", "H2", "H3", "H4"),
                      time = c('t0','t0','t0','t0','t100','t100','t100','t100'),
                      stringsAsFactors = F)
# library(purrr)
df2 <- pmap_df(my_data, ~ data_frame(x = seq(-10, 10, by = 0.1), test = ..3, time = ..4, density = dnorm(x, ..1, ..2)))

ggplot() + 
  geom_ribbon(data=df2[df2$time=='t100',],aes(fill = time, x = x, ymax = density),ymin=0,color='black') +
  geom_ribbon(data=df2[df2$time=='t0',], aes(fill = time, x = x, ymin = -density), ymax=0,color='black') +
  geom_hline(yintercept = 0, color='black')+
  labs(x='fitness')+
  scale_fill_manual(values=c('t0'='grey80','t100'='grey40',labels=c('t0'='t = 0', 't100'='t = 100'),title=''))+
  coord_flip(xlim = c(-5,5),ylim = c(-1,1))+
  facet_grid(.~test)+
  theme_light()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle=90),
        strip.background = element_blank(),
        strip.text = element_blank()) 
