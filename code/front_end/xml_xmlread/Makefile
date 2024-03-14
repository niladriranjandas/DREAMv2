CC		= /data2/nmr/gcc/bin/g++
CFLAGS		= -c -Wall
LDFLAGS		= -L/data2/nmr/gcc/lib -Wl,-rpath,/data2/nmr/gcc/lib64 -lxerces-c
SOURCES		= getFromXml.cpp XmlDomDocument.cpp
INCLUDES	= 
OBJECTS		= $(SOURCES:.cpp=.o)
TARGET		= xmldom

all: $(SOURCES) $(TARGET)

$(TARGET): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -rf $(OBJECTS) $(TARGET)
