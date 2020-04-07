using JSON
using Requests

function downloadModelZipFile(filename, url=""; fileToBeRenamed="", location=pwd())

    # modelName
    modelName = "model.zip"

    # create a temporary directory
    tempDirPath = mktempdir()

    # change to the temporary directory
    cd(tempDirPath)

    println(" > downloadModelZipFile | The directory is: $tempDirPath")

    modelFile = Requests.get(url)
    save(modelFile, modelName)

    # unzip the file
    run(`unzip -qq $modelName`)

    # rename unzipped file if necessary
    if ~isempty(fileToBeRenamed)
        println(" > downloadModelZipFile | Moving $fileToBeRenamed to $filename")
        mv(fileToBeRenamed, filename, remove_destination=true);
    end

    # remove the temporary folder
    rm(tempDirPath, recursive=true)

    # print feedback message
    println(" > downloadModelZipFile | Temporary download folder $tempDirPath removed.")
    println(" > downloadModelZipFile | Model from .zip file stored in $location")
end

# create a temporary directory
tempDirPath = mktempdir()
repoDir = "/var/lib/jenkins/COBRA.models"

# change directory
cd(repoDir)

# parse the JSON list of models
listModels = JSON.parsefile("models.json")

println(" > The temporary directory is: $tempDirPath")

# counters
counter = 0

# loop through the models and download them
for model in listModels
    modelExt = model["name"][end-3:end]
    modelName = tempDirPath*"/"*model["name"]

    if model["url"][end-3:end] == ".mat" || model["url"][end-3:end]  == ".xml"
        # download the file name
        modelFile = get(model["url"])
        save(modelFile, modelName)
        # get the original model
        modelOrig = modelExt[2:end]*"/"*model["name"]

    else
        modelOrig = model["name"][end-2:end]*"/"*model["name"]
        downloadModelZipFile(modelName, model["url"], fileToBeRenamed=model["fileToBeRenamed"], location=tempDirPath)
    end

    # save the checksum
    checkSum = Base.crc32c(read(modelName))
    checkSumOrig = Base.crc32c(read(repoDir*"/"*modelOrig))

    if checkSum == checkSumOrig
        println(" - SAME: $(model["name"]) [new: $checkSum | orig: $checkSumOrig]")
    else
        println(" + DIFFERENT: $(model["name"]) [new: $checkSum | orig: $checkSumOrig]")
        counter = counter + 1
    end
end

# print a feedback message
if counter > 0
    run(pipeline(`echo "[COBRA.models] Models are different"`, `mail -s "+ [COBRA.models] Models are different" laurent.heirendt@uni.lu`))
else
    run(pipeline(`echo "[COBRA.models] Models are the same"`, `mail -s "- [COBRA.models] Models are the same" laurent.heirendt@uni.lu`))
end

# remove the temporary folder
rm(tempDirPath, recursive=true)

# print feedback message
println(" > Folder $tempDirPath has been removed.\n > Done.")

# exit
exit()
