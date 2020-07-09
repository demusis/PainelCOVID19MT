install.packages('rsconnect')
library(rsconnect)

rsconnect::setAccountInfo(name='sesp', 
                          token='8DE82C78F45838E60022AA0CD92E2510', 
                          secret='zY95pmpz5SvJm4szin3XmX5uYrklatR0+YLACTZh')


rsconnect::deployApp('D:/OneDrive/Shiny')

rsconnect::showLogs()
