function cstime = datenum2ticks(t)
cstime = 10^7*60*60*24*(t - 367);