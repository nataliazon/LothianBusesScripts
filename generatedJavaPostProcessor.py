__author__ = 'natalia'

##########
# README
# Look below for the main() function to include the path to your java file.
#########


def addImports(fileLines):
    lastImportLine = 0
    counter = 0
    for line in fileLines:
        if line.find("import")!= -1:
            lastImportLine = counter
        counter +=1
    fileLines.insert(lastImportLine+1, "import java.text.SimpleDateFormat; // Natalia \n")
    fileLines.insert(lastImportLine+2, "import java.util.ArrayList;\n")
    fileLines.insert(lastImportLine+3, "import java.util.Date;\n")
    counter = 0
    for line in fileLines:
        counter +=1
    return fileLines

def addFields(fileLines):
    classStartLine = 0
    counter = 0
    for line in fileLines:
        if line.find("extends CarmaModel")!= -1:
            classStartLine = counter
        counter +=1
    fileLines.insert(classStartLine + 1, "    String timeStamp = new SimpleDateFormat(\"yyyy.MM.dd.HH.mm\").format(new Date());\n")
    fileLines.insert(classStartLine + 2, "    FastAppender fastAppenderMoveUpdate = new FastAppender(\"/Users/natalia/Desktop/\"+timeStamp+\"_moveUpdateD.txt\");\n")
    fileLines.insert(classStartLine + 3, "    FastAppender fastAppenderDepartureTimes = new FastAppender(\"/Users/natalia/Desktop/\"+timeStamp+\"_depC.txt\");\n")

    counter = 0
    for line in fileLines:
        counter +=1

    return fileLines

def addContentToBroadcastUpdateMove(fileLines):

    broadcastUpdateMoveIfClause = 0
    counter = 0
    for line in fileLines:
        if line.find("if (action==__ACT__move) {")!= -1:
            broadcastUpdateMoveIfClause = counter
        counter +=1
    fileLines.insert(broadcastUpdateMoveIfClause + 1, "                Integer __SENDER__myId = (Integer) sender.get( \"myId\" );\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 2, "                List<String> listOfValues = new ArrayList<String>();\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 3, "                \n")
    fileLines.insert(broadcastUpdateMoveIfClause + 4, "                listOfValues.add(\"bId: \" + __SENDER__myId);\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 5, "                listOfValues.add(\"curNow: \" + now);\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 6, "                listOfValues.add(\"bStartTime: \" + __SENDER__startTime);\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 7, "                listOfValues.add(\"bCurLoc: \" + __SENDER__currentLocation);\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 8, "                listOfValues.add(\"bCurLoc_x: \" + get(__CONST__x,__SENDER__currentLocation));\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 9, "                listOfValues.add(\"bCurLoc_y: \" + get(__CONST__y,__SENDER__currentLocation));\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 10, "                \n")
    fileLines.insert(broadcastUpdateMoveIfClause + 11, "                fastAppenderMoveUpdate.appendToFile(listOfValues);\n")
    fileLines.insert(broadcastUpdateMoveIfClause + 12, "                \n")
    return fileLines


def main():

    print("Running post processing ...")

    ### The variable below should hold the path to the java file
    pathToFileForReading = "/Users/natalia/Work/workspace/LothianBuses/src-gen/ms/lothian.java"


    openedFile = open(pathToFileForReading)
    readLines = openedFile.readlines()
    openedFile.close()

    readLines = addImports(readLines)
    readLines = addFields(readLines)
    readLines = addContentToBroadcastUpdateMove(readLines)

    counter = 0
    for line in readLines:
        print(str(counter) + ": " + line)
        counter +=1

    pathToFileForWriting = pathToFileForReading
    openedForWriting = open(pathToFileForWriting, "w")

    for line in readLines:
        openedForWriting.write(line)
    openedForWriting.close()

if __name__ == "__main__":
    main()