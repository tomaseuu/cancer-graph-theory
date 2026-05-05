import os

from dotenv import load_dotenv
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base, sessionmaker


load_dotenv()

DATABASE_URL = os.getenv("DATABASE_URL", "").strip()
DATABASE_ENABLED = bool(DATABASE_URL)

engine = create_engine(DATABASE_URL, pool_pre_ping=True) if DATABASE_ENABLED else None
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine) if DATABASE_ENABLED else None
Base = declarative_base()


def get_db():
    if not DATABASE_ENABLED or SessionLocal is None:
        return

    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db():
    if not DATABASE_ENABLED or engine is None:
        return

    from app.db import models  # noqa: F401

    Base.metadata.create_all(bind=engine)
