from fastapi import FastAPI

from app.api.routes import router
from app.db.database import DATABASE_ENABLED, init_db


app = FastAPI(
    title="OncoGraph API",
    description="Graph-based cancer gene discovery API",
    version="0.1.0",
)


@app.on_event("startup")
def startup_event():
    if DATABASE_ENABLED:
        init_db()


app.include_router(router)
