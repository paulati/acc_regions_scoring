base_path <- "/paula/2018/acc_regions/scoring/data/results/filter/join"
setwd(base_path)

file_name <- "join_filtered_elements.csv"

data <- read.table(file_name, sep="\t", stringsAsFactors = FALSE)
colnames(data) <- c("id", "len", "coord", "shift_count", "hamming_in_out")

#data.max <- data[data$hamming_in_out > 0, ]

data.max <- data

data.max$shift_count_rel <- data.max$shift_count / data.max$len

# escalo a 0 1
max.distance <- max(data.max$hamming_in_out)
min.distance <- min(data.max$hamming_in_out)
data.max$hamming_in_out_scl <- (data.max$hamming_in_out - rep(min.distance, nrow(data.max)))/ (max.distance - min.distance)

max.shift.count <- max(data.max$shift_count_rel)
min.shift.count <- min(data.max$shift_count_rel)
data.max$shift_count_rel_scl <- (data.max$shift_count_rel - rep(min.shift.count, nrow(data.max))) / (max.shift.count - min.shift.count)

plot(data.max$shift_count_rel_scl, data.max$hamming_in_out_scl)
hist(data.max$shift_count_rel_scl)
plot(data.max$shift_count_rel_scl)
hist(data.max$hamming_in_out_scl)

#max.distance <- data.max[(data.max$distance.rel > 0.8), ]
#max.distance.shift <- max.distance[max.distance$shift.rel > 0.6, ]

# calculo el modulo del vector formado por distance.rel y shift.rel
data.max$norm2_scl <- apply(data.max[, c("hamming_in_out_scl", "shift_count_rel_scl")], 1, function(x) norm(x, type="2"))


data.max.sort <- data.max[order(-data.max$norm2_scl), ]

#write.table(data.max.sort, "join_filtered_elements_norm.csv", sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


#hist(data.max$norm2_scl)

# plot(data.max$norm2, data.max$V2)  


