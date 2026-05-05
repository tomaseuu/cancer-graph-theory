from fastapi import FastAPI

from app.api.routes import router


app = FastAPI(
    title="OncoGraph API",
    description="Graph-based cancer gene discovery API",
    version="0.1.0",
)

app.include_router(router)
