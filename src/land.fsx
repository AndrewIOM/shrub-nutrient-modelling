#load "../packages/Bristlecone/bristlecone.fsx"

open FSharp.Data

module Options =
    let saveFolder = "/Users/andrewmartin/Desktop/"

    let east,west,north,south = (1, 1, 1, 1)
    let boxDistance = 256
    let urlBase =
        @"http://visionofbritain.org.uk/geowebcache/service/wms?LAYERS=land&FORMAT=image%2Fjpeg&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&EXCEPTIONS=application%2Fvnd.ogc.se_inimage&SRS=epsg%3A3034"


let bottomCorners =
    let lats = [ Options.east .. Options.boxDistance .. Options.west ]
    let lons = [ Options.north .. Options.boxDistance .. Options.south ]
    lats |> List.zip lons

let coordinateString one two three four =
    sprintf "&BBOX=%i,%i,%i,%i&WIDTH=256&HEIGHT=256" one two three four

let requestTile tileLocation = async {
    let! response = Http.AsyncRequest tileLocation
    if response.StatusCode = 200
    then return Some response.Body
    else return None
}

let requests =
    bottomCorners
    |> List.map(fun (x,y) ->
        let tileLocation = Options.urlBase + (coordinateString x (x + Options.boxDistance) y (y + Options.boxDistance))
        ((x, y), requestTile tileLocation))

for ((x,y),request) in requests do
    let result = request |> Async.RunSynchronously
    match result with
    | Some r -> 
        match r with
        | Binary b -> 
            let fileName = sprintf "cool-tile-%i-%i.jpeg" x y
            System.IO.File.WriteAllBytes(Options.saveFolder + fileName, b)
        | Text s -> invalidOp "Shouldn't be a string, was the request valid??"
    | None -> ()
    System.Threading.Thread.Sleep(10000)
