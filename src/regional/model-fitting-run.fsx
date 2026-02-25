#load "model-fitting.fsx"

``Model-fitting``.work
|> Seq.rev
|> Seq.iter (Bristlecone.Workflow.Orchestration.OrchestrationMessage.StartWorkPackage >> ``Model-fitting``.orchestrator.Post)

System.Console.ReadLine()
